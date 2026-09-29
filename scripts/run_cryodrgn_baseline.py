"""Run the installed cryoDRGN classical voxel backprojection, not neural ab initio.

The software's default real-space window is retained and recorded, so this is
an external sanity check, not the strictly matched main baseline.
"""

import argparse, json, pickle, subprocess, time, os
from pathlib import Path
import numpy as np, mrcfile
from fourier_splats.fsc import fsc, resolution

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("dataset")
p.add_argument("--tag", default="main")
args = p.parse_args()
d = ROOT / "data" / args.dataset
o = ROOT / "results" / args.tag / args.dataset
prep = d / "cryodrgn"
prep.mkdir(exist_ok=True)
meta = np.load(d / "metadata.npz")
manifest = json.loads((d / "manifest.json").read_text())
parts = np.load(o / "partitions.npz")
images = np.load(d / "images.npy", mmap_mode="r")
D = images.shape[-1]
# Match per-image background normalization used by the proposed method.
q = np.arange(-D // 2, D // 2)
y, x = np.meshgrid(q, q, indexing="ij")
outer = x * x + y * y > (D * 0.43) ** 2
scale = np.maximum(np.asarray(images[:, outer]).std(axis=1), 1e-8)
with mrcfile.new(str(prep / "particles.mrcs"), overwrite=True) as m:
    m.set_data((images / scale[:, None, None]).astype("float32"))
    m.voxel_size = manifest["pixel_size_A"]
pickle.dump((meta["rotations"], meta["translations"]), open(prep / "poses.pkl", "wb"))
ct = meta["ctf"].copy()
ct[:, 0] = D
ct[:, 1] = manifest["pixel_size_A"]
pickle.dump(ct, open(prep / "ctf.pkl", "wb"))
commands = []
elapsed = []
volumes = []
for h in [0, 1]:
    pickle.dump(parts[f"half{h}"], open(prep / f"half{h}.pkl", "wb"))
    output = o / f"cryodrgn-backproject-half{h}"
    command = [
        str(ROOT / ".venv/bin/cryodrgn"),
        "backproject_voxel",
        str(prep / "particles.mrcs"),
        "--poses",
        str(prep / "poses.pkl"),
        "--ctf",
        str(prep / "ctf.pkl"),
        "--ind",
        str(prep / f"half{h}.pkl"),
        "-o",
        str(output),
        "--no-half-maps",
        "--ctf-alg",
        "mul",
        "--reg-weight",
        "1.0",
        "-b",
        "64",
    ]
    if manifest["data_sign"] == 1:
        command += ["--uninvert-data"]
    commands.append(command)
    start = time.time()
    with open(ROOT / "logs" / f"cryodrgn-{args.dataset}-half{h}.log", "w") as log:
        subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            env={**os.environ, "OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "4"},
        )
    elapsed.append(time.time() - start)
    with mrcfile.open(output / "backproject.mrc") as m:
        v = m.data.copy()
    volumes.append(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))))
curve = fsc(*volumes, manifest["pixel_size_A"])
np.savetxt(
    o / "cryodrgn-backproject-half-fsc.csv",
    curve,
    delimiter=",",
    header="shell,frequency_inv_A,fsc,voxel_count",
    comments="",
)
gauss = sum(np.load(o / f"gaussian-half{h}.npz")["fourier"] for h in [0, 1]) / 2
cross = fsc(gauss, (volumes[0] + volumes[1]) / 2, manifest["pixel_size_A"])
np.savetxt(
    o / "gaussian-vs-cryodrgn-backproject-fsc.csv",
    cross,
    delimiter=",",
    header="shell,frequency_inv_A,fsc,voxel_count",
    comments="",
)
record = {
    "method": "cryoDRGN 4.3.1 backproject_voxel; classical, fixed-pose",
    "commands": commands,
    "seconds": elapsed,
    "default_real_space_window_retained": True,
    "uses_all_in_band_fourier_pixels": True,
    "half_fsc_resolution": resolution(curve),
    "cross_method_mean_fsc": float(np.nanmean(cross[:, 2])),
    "cross_method_threshold_0.5": resolution(cross, 0.5),
}
(o / "cryodrgn-baseline.json").write_text(json.dumps(record, indent=2))
print(json.dumps(record, indent=2))
