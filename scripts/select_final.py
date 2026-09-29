"""Freeze regularization from RAG validation only, and collect final runs."""

from pathlib import Path
import json, shutil, numpy as np, mrcfile
from fourier_splats.fsc import fsc, resolution

ROOT = Path(__file__).resolve().parents[1]
choice = {}
candidates = {}
for method in ["gaussian", "voxel"]:
    values = {
        str(reg): json.loads(
            (
                ROOT / "results" / ("full-reg" + str(reg)) / "10049/metrics.json"
            ).read_text()
        )["methods"][method]["validation_nmse"]
        for reg in [0.1, 1, 10]
    }
    choice[method] = min(values, key=values.get)
    candidates[method] = values
assert choice == {"gaussian": "1", "voxel": "1"}, choice
selection = {
    "selected_regularization_factor": 1,
    "selection_dataset": "EMPIAR-10049",
    "criterion": "minimum validation image NMSE, 409 held-out particles",
    "candidates": candidates,
    "frozen_for_other_datasets": True,
    "frequency_samples": "all 1410 nonredundant Fourier pixels inside radius 30 at D=64",
    "image_window": "cryoDRGN default linear taper, 0.85 to 0.99 of half-box",
    "important": "The 0.143 FSC did not select this regularizer. Lambda=10 has high FSC but worse validation prediction.",
}
(ROOT / "results/selection.json").write_text(json.dumps(selection, indent=2))
for d in ["10028", "10049", "10076"]:
    src = ROOT / "results/full-reg1" / d
    if not (src / "metrics.json").exists():
        continue
    j = json.loads((src / "metrics.json").read_text())
    if "total_seconds" not in j:
        continue
    dst = ROOT / "results/final" / d
    shutil.copytree(src, dst, dirs_exist_ok=True)
    old = ROOT / "results/main" / d
    if (old / "cryodrgn-baseline.json").exists():
        f = []
        for h in [0, 1]:
            with mrcfile.open(
                old / f"cryodrgn-backproject-half{h}/backproject.mrc"
            ) as m:
                v = m.data.copy()
            F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
            D = v.shape[0]
            q = np.arange(-D // 2, D // 2)
            z, y, x = np.meshgrid(q, q, q, indexing="ij")
            F[x * x + y * y + z * z > (D / 2 - 2) ** 2] = 0
            f.append(F)
            shutil.copy(
                old / f"cryodrgn-backproject-half{h}/backproject.mrc",
                dst / f"cryodrgn-backproject-half{h}.mrc",
            )
        curve = fsc(*f, j["pixel_size_A"])
        np.savetxt(
            dst / "cryodrgn-backproject-half-fsc.csv",
            curve,
            delimiter=",",
            header="shell,frequency_inv_A,fsc,voxel_count",
            comments="",
        )
        gauss = (
            sum(np.load(dst / f"gaussian-half{h}.npz")["fourier"] for h in [0, 1]) / 2
        )
        cross = fsc(gauss, (f[0] + f[1]) / 2, j["pixel_size_A"])
        np.savetxt(
            dst / "gaussian-vs-cryodrgn-backproject-fsc.csv",
            cross,
            delimiter=",",
            header="shell,frequency_inv_A,fsc,voxel_count",
            comments="",
        )
        rec = json.loads((old / "cryodrgn-baseline.json").read_text())
        rec["half_fsc_resolution"] = resolution(curve)
        rec["cross_method_mean_fsc"] = float(cross[:, 2].mean())
        rec["cross_method_threshold_0.5"] = resolution(cross, 0.5)
        rec["analysis_sphere_radius"] = D / 2 - 2
        (dst / "cryodrgn-baseline.json").write_text(json.dumps(rec, indent=2))
    # Paired uncertainty interval, treating particles as sampling units (not micrographs).
    a = np.load(dst / "gaussian-test-per-particle-mse.npy")
    b = np.load(dst / "voxel-test-per-particle-mse.npy")
    diff = a - b
    rng = np.random.default_rng(992)
    bootstrap = diff[rng.integers(len(diff), size=(5000, len(diff)))].mean(axis=1)
    summary = {
        "mean_gaussian_minus_voxel_mse": float(diff.mean()),
        "bootstrap_95_percent_interval": np.quantile(
            bootstrap, [0.025, 0.975]
        ).tolist(),
        "bootstrap_unit": "particle; conditional on supplied poses; micrograph correlations not accounted for",
        "replicates": 5000,
    }
    (dst / "paired-test-bootstrap.json").write_text(json.dumps(summary, indent=2))
    print("FINAL", d, j["methods"]["gaussian"]["half_fsc_resolution"])
