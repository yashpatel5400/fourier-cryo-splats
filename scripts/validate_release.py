"""Validate saved maps, particle partitions and published numerical curves."""

from pathlib import Path
import json, numpy as np, mrcfile
from fourier_splats.physics import volume_from_fourier
from fourier_splats.fsc import fsc

root = Path(__file__).resolve().parents[1]
records = {}
for d in ["10028", "10049", "10076"]:
    p = root / "results/final" / d
    j = json.loads((p / "metrics.json").read_text())
    partition = np.load(p / "partitions.npz")
    ids = np.concatenate(
        [partition[key] for key in ["half0", "half1", "validation", "test"]]
    )
    assert np.array_equal(np.sort(ids), np.arange(8192))
    assert len(np.unique(partition["source_indices"])) == 8192
    for method in ["gaussian", "voxel"]:
        fourier = []
        for h in [0, 1]:
            F = np.load(p / f"{method}-half{h}.npz")["fourier"]
            fourier.append(F)
            with mrcfile.open(p / f"{method}-half{h}.mrc") as m:
                assert m.data.shape == (64, 64, 64)
                np.testing.assert_allclose(
                    m.data, volume_from_fourier(F), rtol=1e-6, atol=1e-8
                )
                np.testing.assert_allclose(
                    float(m.voxel_size.x), j["pixel_size_A"], rtol=1e-6
                )
        expected = fsc(*fourier, j["pixel_size_A"])
        stored = np.loadtxt(p / f"{method}-half-fsc.csv", delimiter=",", skiprows=1)
        np.testing.assert_allclose(stored, expected, rtol=1e-8, atol=1e-8)
        for half in j["methods"][method]["halves"]:
            for component in half.values():
                assert (
                    component["cg_info"] == 0
                    and component["relative_normal_residual"] < 1e-5
                )
    records[d] = {
        "disjoint_partitions": True,
        "unique_source_particles": 8192,
        "map_checkpoint_match": True,
        "physical_pixel_size_verified": True,
        "FSC_recomputed_matches_CSV": True,
        "all_CG_solves_converged": True,
    }
(root / "results/validation/release-checks.json").write_text(
    json.dumps(records, indent=2)
)
print(json.dumps(records, indent=2))
