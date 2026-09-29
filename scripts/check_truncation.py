"""Check Gaussian rendering cutoff against a larger neighborhood at fixed coefficients."""

import json
from pathlib import Path
import numpy as np
from fourier_splats.basis import design

root = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(123)
k = rng.uniform(-30, 30, (60000, 3)).astype("float32")
k = k[np.linalg.norm(k, axis=1) < 30][:20000]
records = {}
for d in ["10028", "10049", "10076"]:
    p = root / "results/final" / d / "gaussian-half0.npz"
    if not p.exists():
        continue
    a = np.load(p)
    pred = []
    for radius in [2.0, 3.0]:
        ar, ai = design(k, 64, "gaussian", 0.5, radius)
        pred.append(ar @ a["real"] + 1j * (ai @ a["imag"]))
    relative = float(np.linalg.norm(pred[0] - pred[1]) / np.linalg.norm(pred[1]))
    records[d] = {
        "relative_complex_L2_error": relative,
        "evaluation_points": len(k),
        "reference_radius": 3,
        "implemented_radius": 2,
        "scope": "fixed-coefficient rendering comparison, not refitting or a resolution measurement",
    }
    assert relative < 0.005, (d, relative)
(root / "results/validation/truncation.json").write_text(json.dumps(records, indent=2))
print(records)
