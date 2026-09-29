#!/usr/bin/env python3
"""Fit true experimental image likelihoods and save all FSC/prediction diagnostics."""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import argparse, csv, json, time, platform, gc
from pathlib import Path
import numpy as np
import mrcfile
from scipy.sparse import vstack
from scipy.sparse.linalg import LinearOperator, cg
from fourier_splats.physics import fft_center, ctf, volume_from_fourier
from fourier_splats.basis import design, evaluate_grid
from fourier_splats.fsc import fsc, resolution

ROOT = Path(__file__).resolve().parents[1]


def observations(directory, box, samples, seed, window=False):
    images = np.load(directory / "images.npy", mmap_mode="r")
    meta = np.load(directory / "metadata.npz")
    manifest = json.loads((directory / "manifest.json").read_text())
    original = images.shape[-1]
    q = np.arange(-box // 2, box // 2, dtype=np.float32)
    y, x = np.meshgrid(q, q, indexing="ij")
    mask = ((x * x + y * y) <= (box / 2 - 2) ** 2) & ((y > 0) | ((y == 0) & (x > 0)))
    coords = np.stack([x[mask], y[mask]], axis=-1)
    if samples > len(coords):
        samples = len(coords)
    rng = np.random.default_rng(seed)
    n = len(images)
    k = np.empty((n, samples, 3), np.float32)
    obs = np.empty((n, samples), np.complex64)
    transfer = np.empty((n, samples), np.float32)
    normalization = np.empty(n, np.float32)
    for start in range(0, n, 32):
        stop = min(n, start + 32)
        raw = np.asarray(images[start:stop])
        if window:
            from cryodrgn.masking import spherical_window_mask

            win = spherical_window_mask(D=original, in_rad=0.85, out_rad=0.99).numpy()
            F = fft_center(raw * win[None, :, :])
        else:
            F = fft_center(raw)
        c = original // 2
        h = box // 2
        F = F[:, c - h : c + h, c - h : c + h]
        # One scale for each image, calculated in the spatial background of the downloaded image.
        qr = np.arange(-original // 2, original // 2)
        yy, xx = np.meshgrid(qr, qr, indexing="ij")
        outer = xx * xx + yy * yy > (original * 0.43) ** 2
        scale = np.maximum(raw[:, outer].std(axis=1), 1e-8)
        normalization[start:stop] = scale
        F = manifest["data_sign"] * F / scale[:, None, None]
        F = F[:, mask]
        shifts = meta["translations"][start:stop]
        # cryoDRGN translates its input images by the stored shift to center them.
        F *= np.exp(-2j * np.pi * np.einsum("pc,bc->bp", coords, shifts))
        C = ctf(
            coords / (manifest["raw_box"] * manifest["raw_pixel_size_A"]),
            meta["ctf"][start:stop],
        )
        selected = np.stack(
            [
                rng.choice(len(coords), samples, replace=False)
                for _ in range(stop - start)
            ]
        )
        qsel = coords[selected]
        plane = np.concatenate(
            [qsel, np.zeros((*qsel.shape[:-1], 1), np.float32)], axis=-1
        )
        k[start:stop] = np.einsum("bpi,bij->bpj", plane, meta["rotations"][start:stop])
        obs[start:stop] = np.take_along_axis(F, selected, axis=1)
        transfer[start:stop] = np.take_along_axis(C, selected, axis=1)
    # Global scaling changes coefficient units only; derived from training images below.
    return k, obs, transfer, manifest, normalization


def matrices(k, c, box, kind, sigma, radius):
    ars = []
    ais = []
    for start in range(0, len(k), 4096):
        a, b = design(k[start : start + 4096], box, kind, sigma, radius)
        w = c[start : start + 4096, None]
        ars.append(a.multiply(w).tocsr())
        ais.append(b.multiply(w).tocsr())
    ar = vstack(ars, format="csr")
    ai = vstack(ais, format="csr")
    return ar, ai


def solve(a, y, reg, iterations):
    t = time.time()
    at = a.T.tocsr()
    diag = np.asarray(a.power(2).sum(axis=0)).ravel()
    positive = diag[diag > 1e-6]
    ridge = reg * np.median(positive)
    d = diag + ridge
    rhs = at @ y
    op = LinearOperator(
        (a.shape[1],) * 2, matvec=lambda x: at @ (a @ x) + ridge * x, dtype=np.float32
    )
    pre = LinearOperator(op.shape, matvec=lambda x: x / d, dtype=np.float32)
    history = []

    def callback(x):
        if len(history) % 10 == 0:
            residual = op @ x - rhs
            history.append(
                float(np.linalg.norm(residual) / max(np.linalg.norm(rhs), 1e-20))
            )
        else:
            history.append(None)

    solution, info = cg(
        op, rhs, rtol=1e-5, atol=0, maxiter=iterations, M=pre, callback=callback
    )
    residual = float(
        np.linalg.norm(op @ solution - rhs) / max(np.linalg.norm(rhs), 1e-20)
    )
    record = {
        "regularization_factor": reg,
        "ridge_lambda": float(ridge),
        "cg_info": int(info),
        "iterations": len(history),
        "relative_normal_residual": residual,
        "seconds": time.time() - t,
        "observations": a.shape[0],
        "unknowns": a.shape[1],
        "nonzeros": a.nnz,
        "residual_history": history,
    }
    return solution, record


def predict(k, c, real, imag, args, kind):
    pred = np.empty(len(k), np.complex64)
    for start in range(0, len(k), 8192):
        ar, ai = design(
            k[start : start + 8192], args.box, kind, args.sigma, args.radius
        )
        pred[start : start + len(ar.indptr) - 1] = (ar @ real + 1j * (ai @ imag)) * c[
            start : start + len(ar.indptr) - 1
        ]
    return pred


def save_map(path, f, pixel_size):
    with mrcfile.new(str(path), overwrite=True) as m:
        m.set_data(volume_from_fourier(f))
        m.voxel_size = pixel_size


def run(args):
    directory = ROOT / "data" / args.dataset
    out = ROOT / "results" / args.tag / args.dataset
    out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    k, y, c, manifest, normalization = observations(
        directory, args.box, args.samples, args.seed, args.window
    )
    if args.phase_randomize:
        y *= np.exp(
            2j * np.pi * np.random.default_rng(args.seed + 999).uniform(size=y.shape)
        )
    n = len(y)
    rng = np.random.default_rng(args.seed + 1)
    order = rng.permutation(n)
    ntest = max(1, n // 20)
    test = order[:ntest]
    validation = order[ntest : 2 * ntest]
    train = order[2 * ntest :]
    halves = [train[::2], train[1::2]]
    if args.group_splits:
        rows=list(csv.DictReader((ROOT/'research/uncertainty/splits'/(args.dataset+'.csv')).open()))
        by_split={label:np.array([int(r['output_index']) for r in rows if r['split']==label])
                  for label in ['inference_half0','inference_half1','tune','test']}
        halves=[by_split['inference_half0'],by_split['inference_half1']]
        train=np.concatenate(halves);validation=by_split['tune'];test=by_split['test']
    scale = np.sqrt(np.mean(np.abs(y[train]) ** 2))
    y /= scale
    np.savez(
        out / "partitions.npz",
        half0=halves[0],
        half1=halves[1],
        validation=validation,
        test=test,
        source_indices=np.load(directory / "indices.npy"),
    )
    print(
        "OBSERVATIONS",
        args.dataset,
        "box",
        args.box,
        "train",
        len(train),
        "samples/image",
        y.shape[1],
        flush=True,
    )
    metrics = {
        "dataset": args.dataset,
        "particle_count": n,
        "train_particles": len(train),
        "validation_particles": len(validation),
        "test_particles": len(test),
        "half_particles": [len(x) for x in halves],
        "box": args.box,
        "pixel_size_A": manifest["raw_box"] * manifest["raw_pixel_size_A"] / args.box,
        "samples_per_image": y.shape[1],
        "seed": args.seed,
        "sigma": args.sigma,
        "radius": args.radius,
        "data_sign": manifest["data_sign"],
        "window": args.window,
        "phase_randomize": args.phase_randomize,
        "source_group_splits": args.group_splits,
        "regularization_factor": args.reg,
        "pose_protocol": "fixed public consensus poses; conditional half-map FSC, not fully gold-standard",
        "methods": {},
        "platform": platform.platform(),
    }
    allmaps = {}
    for kind in args.methods.split(","):
        if kind not in ["gaussian", "voxel"]:
            raise ValueError(kind)
        tm = time.time()
        method = {"halves": []}
        maps = []
        coeff = []
        for h, ids in enumerate(halves):
            print("BUILD", kind, "half", h, flush=True)
            tb = time.time()
            ar, ai = matrices(
                k[ids].reshape(-1, 3),
                c[ids].ravel(),
                args.box,
                kind,
                args.sigma,
                args.radius,
            )
            print(
                "SOLVE",
                kind,
                "half",
                h,
                "nnz",
                ar.nnz,
                "build seconds",
                round(time.time() - tb, 2),
                flush=True,
            )
            re, lr = solve(ar, y[ids].real.ravel(), args.reg, args.iterations)
            im, li = solve(ai, y[ids].imag.ravel(), args.reg, args.iterations)
            del ar, ai
            gc.collect()
            F = evaluate_grid(re, im, args.box, kind, args.sigma, args.radius)
            maps.append(F)
            coeff.append((re, im))
            save_map(out / f"{kind}-half{h}.mrc", F, metrics["pixel_size_A"])
            np.savez(out / f"{kind}-half{h}.npz", real=re, imag=im, fourier=F)
            method["halves"].append({"real": lr, "imag": li})
            print(
                "DONE",
                kind,
                "half",
                h,
                "residual",
                lr["relative_normal_residual"],
                li["relative_normal_residual"],
                flush=True,
            )
        mean = (maps[0] + maps[1]) / 2
        allmaps[kind] = mean
        curve = fsc(maps[0], maps[1], metrics["pixel_size_A"])
        np.savetxt(
            out / f"{kind}-half-fsc.csv",
            curve,
            delimiter=",",
            header="shell,frequency_inv_A,fsc,voxel_count",
            comments="",
        )
        method["half_fsc_resolution"] = resolution(curve)
        re = (coeff[0][0] + coeff[1][0]) / 2
        im = (coeff[0][1] + coeff[1][1]) / 2
        for label, ids in [("validation", validation), ("test", test)]:
            prediction = predict(
                k[ids].reshape(-1, 3), c[ids].ravel(), re, im, args, kind
            )
            truth = y[ids].ravel()
            per_particle = np.mean(
                np.abs(prediction.reshape(len(ids), -1) - y[ids]) ** 2, axis=1
            )
            method[label + "_nmse"] = float(
                np.mean(np.abs(prediction - truth) ** 2) / np.mean(np.abs(truth) ** 2)
            )
            method[label + "_correlation"] = float(
                np.real(np.vdot(prediction, truth))
                / np.sqrt(
                    np.vdot(prediction, prediction).real * np.vdot(truth, truth).real
                )
            )
            np.save(out / f"{kind}-{label}-per-particle-mse.npy", per_particle)
        method["seconds"] = time.time() - tm
        method["estimated_coefficient_parameters"] = 2 * len(re) - 1
        metrics["methods"][kind] = method
        save_map(out / f"{kind}-mean.mrc", mean, metrics["pixel_size_A"])
        (out / "metrics.json").write_text(json.dumps(metrics, indent=2))
        print(
            "METRICS",
            kind,
            json.dumps({kk: vv for kk, vv in method.items() if kk != "halves"}),
            flush=True,
        )
    if "gaussian" in allmaps and "voxel" in allmaps:
        cross = fsc(allmaps["gaussian"], allmaps["voxel"], metrics["pixel_size_A"])
        np.savetxt(
            out / "gaussian-vs-voxel-fsc.csv",
            cross,
            delimiter=",",
            header="shell,frequency_inv_A,fsc,voxel_count",
            comments="",
        )
        metrics["method_agreement_fsc_threshold_0.5"] = resolution(cross, 0.5)
        metrics["mean_method_agreement_fsc"] = float(np.nanmean(cross[:, 2]))
    metrics["total_seconds"] = time.time() - t
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2))
    print("COMPLETE", args.dataset, round(time.time() - t, 1), "seconds", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("dataset")
    p.add_argument("--box", type=int, default=64)
    p.add_argument("--samples", type=int, default=256)
    p.add_argument("--sigma", type=float, default=0.5)
    p.add_argument("--radius", type=float, default=2.0)
    p.add_argument("--reg", type=float, default=0.1)
    p.add_argument("--iterations", type=int, default=80)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--tag", default="main")
    p.add_argument("--methods", default="gaussian,voxel")
    p.add_argument("--window", action="store_true")
    p.add_argument("--phase-randomize", action="store_true")
    p.add_argument("--group-splits", action="store_true",help="Use the audited exposure-group uncertainty development splits")
    run(p.parse_args())
