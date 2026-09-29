"""Verify that download recovery preserves particle-to-metadata identity exactly."""

from pathlib import Path
import hashlib, json, numpy as np

root = Path(__file__).resolve().parents[1]
p = root / "results/validation/source-identity-before-recovery.json"
if p.exists():
    expected = json.loads(p.read_text())
    images = np.load(root / "data/10028/images.npy", mmap_mode="r")
    indices = np.load(root / "data/10028/indices.npy")
    lookup = {str(int(source)): i for i, source in enumerate(indices)}
    count = 0
    for source, digest in expected.items():
        assert source in lookup, (
            "Previously verified source lost during recovery",
            source,
        )
        got = hashlib.sha256(images[lookup[source]].tobytes()).hexdigest()
        assert got == digest, ("Particle-to-metadata identity changed", source)
        count += 1
    result = {
        "verified_preserved_source_identities": count,
        "sha256_per_image_matches": True,
        "selection_count": len(indices),
        "new_selection_indices_unique": len(np.unique(indices)) == len(indices),
    }
    (root / "results/validation/source-recovery-check.json").write_text(
        json.dumps(result, indent=2)
    )
    print(result)
else:
    print("Fresh download: no interrupted-run recovery audit is applicable.")
