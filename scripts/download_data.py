#!/usr/bin/env python3
"""Fetch actual deposited particle bytes with reproducible, verified HTTP ranges.

The selection is stratified contiguous blocks over the original particle ordering,
then the published cryoDRGN particle filter (where supplied) is applied. This is
not a uniform independent sample and is explicitly recorded in the manifest.
"""

import argparse, concurrent.futures as cf, hashlib, json, pickle, struct, time
import os
from pathlib import Path
import numpy as np
import requests
from scipy.fft import fft2, fftshift, ifftshift, ifft2

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://ftp.ebi.ac.uk/empiar/world_availability"


def get_range(url, start, end):
    for attempt in range(6):
        try:
            # A range-specific query avoids intermediate caches reusing a header-only request.
            response = requests.get(
                url + f"?range={start}-{end}",
                headers={"Range": f"bytes={start}-{end}"},
                timeout=(20, 180),
            )
            response.raise_for_status()
            if (
                response.status_code != 206
                or response.headers.get("Content-Range", "").split("/")[0]
                != f"bytes {start}-{end}"
            ):
                raise RuntimeError(
                    f"Range not honored: {response.status_code}, {response.headers.get('Content-Range')}"
                )
            if len(response.content) != end - start + 1:
                raise RuntimeError("Truncated download")
            return (
                response.content,
                response.headers.get("ETag"),
                response.headers.get("Content-Range"),
            )
        except Exception:
            if attempt == 5:
                raise
            time.sleep(2**attempt)


def run(args):
    directory = Path(args.output_directory).resolve() if args.output_directory else ROOT / "data" / args.dataset
    directory.mkdir(parents=True, exist_ok=True)
    source = ROOT / "background/cryodrgn_empiar" / ("empiar" + args.dataset) / "inputs"
    cs = np.load(next(source.glob("*.cs")), allow_pickle=False)
    # Author-provided metadata repository, never arbitrary user-supplied pickle files.
    poses = pickle.load(open(source / "poses.pkl", "rb"))
    ctf = pickle.load(open(source / "ctf.pkl", "rb"))
    n = len(cs)
    rng = np.random.default_rng(args.seed)
    block = 64
    starts = np.arange(0, n, block)
    rng.shuffle(starts)
    valid = np.ones(n, dtype=bool)
    if (source / "filtered.ind.pkl").exists():
        valid[:] = False
        valid[
            np.asarray(pickle.load(open(source / "filtered.ind.pkl", "rb")), dtype=int)
        ] = True
    published_count = int(valid.sum())
    exclusions_path = ROOT / "provenance" / f"{args.dataset}-source-exclusions.json"
    exclusions = (
        json.loads(exclusions_path.read_text()) if exclusions_path.exists() else []
    )
    for exclusion in exclusions:
        invalid = np.array(
            [
                path.decode().endswith(exclusion["metadata_path_suffix"])
                for path in cs["blob/path"]
            ]
        )
        valid[invalid] = False
    selected = []
    for start in starts:
        selected.extend([i for i in range(start, min(start + block, n)) if valid[i]])
        if len(selected) >= args.count:
            break
    selected = np.sort(selected[: args.count])
    if args.indices_file:
        selected = np.asarray(np.load(args.indices_file, allow_pickle=False), dtype=np.int64)
        if (selected.ndim != 1 or len(selected) != args.count or len(np.unique(selected)) != len(selected)
                or np.any(selected < 0) or np.any(selected >= n) or not np.all(valid[selected])):
            raise ValueError("Explicit selection fails count/uniqueness/source eligibility checks")
        if np.any(np.diff(selected) <= 0):
            raise ValueError("Explicit source indices must be strictly increasing")
    progress_file = directory / "download-progress.json"
    records = []
    if args.resume and progress_file.exists():
        old_indices = np.load(directory / "indices.npy")
        records = json.loads(progress_file.read_text())
        if not np.array_equal(old_indices, selected):
            if not exclusions:
                raise ValueError(
                    "Resume selection differs without a documented source exclusion"
                )
            old_images = np.load(directory / "images.npy", mmap_mode="r")
            new_images = np.lib.format.open_memmap(
                directory / "images-reindexed.npy",
                mode="w+",
                dtype="float32",
                shape=(len(selected), args.box, args.box),
            )
            new_positions = {
                int(source_index): row for row, source_index in enumerate(selected)
            }
            remapped = []
            for record in records:
                old_rows = record["output_indices"]
                source_rows = old_indices[old_rows]
                if not all(int(i) in new_positions for i in source_rows):
                    continue
                new_rows = [new_positions[int(i)] for i in source_rows]
                new_images[new_rows] = old_images[old_rows]
                record["output_indices"] = new_rows
                remapped.append(record)
            new_images.flush()
            del old_images, new_images
            os.replace(directory / "images-reindexed.npy", directory / "images.npy")
            records = remapped
    np.save(directory / "indices.npy", selected)
    N = len(selected)
    D = args.box
    images = np.lib.format.open_memmap(
        directory / "images.npy",
        mode="r+" if records else "w+",
        dtype="float32",
        shape=(N, D, D),
    )
    if images.shape != (N, D, D):
        raise ValueError("Resume image dimensions differ from the requested dimensions")
    np.savez(
        directory / "metadata.npz",
        rotations=poses[0][selected],
        translations=poses[1][selected],
        ctf=ctf[selected],
        indices=selected,
    )
    jobs = []
    for path in np.unique(cs["blob/path"][selected]):
        positions = np.flatnonzero(cs["blob/path"][selected] == path)
        idx = cs["blob/idx"][selected[positions]].astype(int)
        order = np.argsort(idx)
        idx = idx[order]
        positions = positions[order]
        relative = path.decode().split("/imported/")[-1]
        if args.dataset == "10028":
            relative = "Particles/" + relative
        url = f"{BASE}/{args.dataset}/data/{relative}"
        splits = np.flatnonzero(np.diff(idx) > 16) + 1
        for segment in np.split(np.arange(len(idx)), splits):
            # Bound each response to at most ~64MB for low peak memory.
            for part in np.array_split(
                segment, max(1, int(np.ceil(len(segment) / 128)))
            ):
                if len(part):
                    jobs.append((url, idx[part], positions[part]))
    headers = {}

    def work(job):
        url, ids, pos = job
        if url not in headers:
            raw, _, _ = get_range(url, 0, 1023)
            headers[url] = raw
        header = headers[url]
        nx, ny, nz, mode = struct.unpack("<4i", header[:16])
        nsymbt = struct.unpack("<i", header[92:96])[0]
        if mode != 2 or nx != ny:
            raise ValueError((nx, ny, nz, mode))
        if ids.max() >= nz:
            raise ValueError("Particle index exceeds source stack")
        stride = nx * ny * 4
        start = 1024 + nsymbt + int(ids.min()) * stride
        end = 1024 + nsymbt + (int(ids.max()) + 1) * stride - 1
        raw, etag, content_range = get_range(url, start, end)
        x = np.frombuffer(raw, dtype="<f4").reshape(-1, ny, nx)[ids - ids.min()]
        # Crop the centered Fourier transform with symmetric Nyquist handling.
        f = fftshift(fft2(ifftshift(x, axes=(-2, -1)), workers=1), axes=(-2, -1))
        c = nx // 2
        h = D // 2
        f = f[:, c - h : c + h, c - h : c + h].copy()
        # Remove Nyquist row/column, which otherwise require folding both endpoints.
        f[:, 0, :] = 0
        f[:, :, 0] = 0
        xsmall = (
            fftshift(ifft2(ifftshift(f, axes=(-2, -1)), workers=1), axes=(-2, -1)).real
            * (D / nx) ** 2
        )
        images[pos] = xsmall.astype("float32")
        return {
            "url": url,
            "start": start,
            "end": end,
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "etag": etag,
            "content_range": content_range,
            "header_sha256": hashlib.sha256(header).hexdigest(),
            "source_shape": [nz, ny, nx],
            "source_indices": ids.tolist(),
            "output_indices": pos.tolist(),
        }

    total_jobs = len(jobs)
    job_keys = {(j[0], tuple(j[1].tolist())) for j in jobs}
    records = [r for r in records if (r["url"], tuple(r["source_indices"])) in job_keys]
    completed = {(r["url"], tuple(r["source_indices"])) for r in records}
    jobs = [j for j in jobs if (j[0], tuple(j[1].tolist())) not in completed]
    failures = []
    t = time.time()
    print(
        args.dataset,
        "resuming",
        len(records),
        "verified ranges;",
        len(jobs),
        "remaining",
        flush=True,
    )
    with cf.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(work, j): j for j in jobs}
        for future in cf.as_completed(futures):
            try:
                records.append(future.result())
            except Exception as error:
                job = futures[future]
                failure = {
                    "url": job[0],
                    "indices": job[1].tolist(),
                    "error": repr(error),
                }
                failures.append(failure)
                print("FAILED RANGE", failure, flush=True)
                (directory / "download-failures.json").write_text(
                    json.dumps(failures, indent=2)
                )
                continue
            images.flush()
            if len(records) % 10 == 0 or len(records) == total_jobs:
                print(
                    args.dataset,
                    len(records),
                    "/",
                    total_jobs,
                    "ranges;",
                    round(sum(x["bytes"] for x in records) / 1e9, 2),
                    "GB;",
                    round(time.time() - t),
                    "s",
                    flush=True,
                )
            (directory / "download-progress.json").write_text(json.dumps(records))
    if failures:
        raise RuntimeError(
            f"{len(failures)} download ranges failed; verified progress saved. Retry with --resume."
        )
    manifest = {
        "dataset": args.dataset,
        "count": N,
        "available_particles": n,
        "published_filtered_count": published_count,
        "eligible_after_source_exclusions": int(valid.sum()),
        "box": D,
        "seed": args.seed,
        "selection": "explicit frozen indices" if args.indices_file else "randomly permuted contiguous source blocks of 64, applying published filter; sorted before output",
        "selection_file_sha256": hashlib.sha256(Path(args.indices_file).read_bytes()).hexdigest() if args.indices_file else None,
        "raw_pixel_size_A": float(cs["blob/psize_A"][0]),
        "raw_box": int(cs["blob/shape"][0, 0]),
        "pixel_size_A": float(cs["blob/psize_A"][0]) * int(cs["blob/shape"][0, 0]) / D,
        "data_sign": 1 if args.dataset == "10076" else -1,
        "source": "experimental extracted particle images; not raw movies",
        "metadata_repository": "https://github.com/zhonge/cryodrgn_empiar",
        "source_exclusions": exclusions,
        "elapsed_seconds": time.time() - t,
        "resumed_from_ranges": len(completed),
        "downloaded_bytes": sum(x["bytes"] for x in records),
        "payload_accounting": "verified unique byte ranges; interrupted attempts, retries, HTTP headers and duplicate recovery transfers are not included",
        "ranges": records,
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print("COMPLETE", args.dataset, N, "particles", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("dataset", choices=["10028", "10049", "10076"])
    p.add_argument("--count", type=int, default=8192)
    p.add_argument("--box", type=int, default=64)
    p.add_argument("--seed", type=int, default=20260928)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--indices-file", help="Prespecified eligible source indices in a numeric NPY array")
    p.add_argument("--output-directory", help="Separate output pool; default is data/DATASET")
    run(p.parse_args())
