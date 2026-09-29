import argparse
from pathlib import Path
import numpy as np, matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mrcfile

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--tag", default="main")
args = p.parse_args()
datasets = [
    d
    for d in ["10028", "10049", "10076"]
    if (ROOT / "results" / args.tag / d / "metrics.json").exists()
]
colors = {"gaussian": "#167c80", "voxel": "#d47432", "cryodrgn-backproject": "#5667aa"}
fig, axs = plt.subplots(
    1, len(datasets), figsize=(4.2 * len(datasets), 3.4), squeeze=False
)
for ax, d in zip(axs[0], datasets):
    path = ROOT / "results" / args.tag / d
    for name, label in [
        ("gaussian", "Fourier Gaussian"),
        ("voxel", "Matched voxel"),
        ("cryodrgn-backproject", "cryoDRGN backprojection"),
    ]:
        file = path / (name + "-half-fsc.csv")
        if file.exists():
            a = np.loadtxt(file, delimiter=",", skiprows=1)
            ax.plot(a[:, 1], a[:, 2], label=label, color=colors[name], lw=1.8)
    ax.axhline(0.143, color="grey", ls="--", lw=0.8)
    ax.set(
        title="EMPIAR-" + d,
        xlabel="Spatial frequency (1/Å)",
        ylabel="Conditional half-map FSC",
        ylim=(-0.1, 1.04),
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        frameon=True, framealpha=0.95, edgecolor="none", fontsize=7.5, loc="lower left"
    )
fig.tight_layout()
(ROOT / "figures").mkdir(exist_ok=True)
fig.savefig(ROOT / "figures" / f"{args.tag}-fsc.png", dpi=180)
plt.close(fig)
fig, axs = plt.subplots(
    len(datasets), 6, figsize=(12, 2.15 * len(datasets)), squeeze=False
)
for row, d in enumerate(datasets):
    path = ROOT / "results" / args.tag / d
    for methodindex, method in enumerate(["gaussian", "voxel"]):
        with mrcfile.open(path / (method + "-mean.mrc")) as m:
            v = m.data.copy()
        # Low-pass for illustration only; FSC always uses unfiltered Fourier maps.
        from scipy.ndimage import gaussian_filter

        v = gaussian_filter(v, 1)
        for axis in range(3):
            image = np.take(v, v.shape[axis] // 2, axis=axis)
            ax = axs[row, methodindex * 3 + axis]
            lims = np.quantile(v, [0.01, 0.998])
            ax.imshow(image, cmap="magma", vmin=lims[0], vmax=lims[1], origin="lower")
            ax.set_xticks([])
            ax.set_yticks([])
            if row == 0:
                ax.set_title(
                    method.title() + " " + ["xy", "xz", "yz"][axis], fontsize=9
                )
        axs[row, 0].set_ylabel("EMPIAR-" + d)
fig.tight_layout()
fig.savefig(ROOT / "figures" / f"{args.tag}-slices.png", dpi=180)
plt.close(fig)
