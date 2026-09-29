from pathlib import Path
import json, hashlib, subprocess, numpy as np

root = Path(__file__).resolve().parents[1]
out = root / "provenance"
out.mkdir(exist_ok=True)
source_commit = subprocess.check_output(
    ["git", "-C", str(root / "background/cryodrgn_empiar"), "rev-parse", "HEAD"],
    text=True,
).strip()
for d in ["10028", "10049", "10076"]:
    path = root / "data" / d
    if not (path / "manifest.json").exists():
        continue
    manifest = json.loads((path / "manifest.json").read_text())
    manifest["metadata_commit"] = source_commit
    for name in ["images.npy", "metadata.npz", "indices.npy"]:
        h = hashlib.sha256()
        with open(path / name, "rb") as f:
            for b in iter(lambda: f.read(4 * 1024 * 1024), b""):
                h.update(b)
        manifest[name + "_sha256"] = h.hexdigest()
    images = np.load(path / "images.npy", mmap_mode="r")
    assert np.isfinite(images).all() and (images.std(axis=(1, 2)) > 0).all()
    indices = np.load(path / "indices.npy")
    assert len(np.unique(indices)) == len(indices)
    assert sum(len(x["output_indices"]) for x in manifest["ranges"]) == len(indices)
    out.joinpath(d + "-manifest.json").write_text(json.dumps(manifest, indent=2))
    np.savetxt(
        out / (d + "-source-indices.csv"),
        indices,
        fmt="%d",
        header="zero_based_source_particle_index",
        comments="",
    )
    print(
        d,
        manifest["count"],
        manifest["downloaded_bytes"],
        "bytes, verified finite and unique",
    )
