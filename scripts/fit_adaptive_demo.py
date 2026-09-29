"""Small real-data nonlinear demonstration; explicitly separate from the main study.

Fits coefficients, centers and full anisotropic precisions by stochastic gradient
on low-frequency observations from one real dataset. No FSC claim is based on it.
"""

import argparse, sys, time, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, torch
from run_experiment import observations
from fourier_splats.adaptive import FourierGaussianPairs

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--dataset", default="10049")
p.add_argument("--steps", type=int, default=1500)
args = p.parse_args()
torch.manual_seed(82)
torch.set_num_threads(4)
rng = np.random.default_rng(82)
k, y, c, manifest, _ = observations(ROOT / "data" / args.dataset, 16, 128, 82, True)
ids = rng.permutation(len(k))
train = ids[:-512]
val = ids[-512:]
scale = np.sqrt(np.mean(np.abs(y[train]) ** 2))
y /= scale
# One member per integer-frequency conjugate pair, to r=7.
z, yy, x = np.meshgrid(
    np.arange(-7, 8), np.arange(-7, 8), np.arange(-7, 8), indexing="ij"
)
centers = np.stack([x.ravel(), yy.ravel(), z.ravel()], -1).astype("float32")
choose = (
    (centers[:, 2] > 0)
    | ((centers[:, 2] == 0) & (centers[:, 1] > 0))
    | ((centers[:, 2] == 0) & (centers[:, 1] == 0) & (centers[:, 0] > 0))
)
centers = centers[choose & (np.linalg.norm(centers, axis=1) <= 7)]
coeff = (
    rng.normal(size=len(centers)) * 1e-3 + 1j * rng.normal(size=len(centers)) * 1e-3
).astype("complex64")
L = np.tile(np.eye(3, dtype="float32") / 0.7, (len(centers), 1, 1))
device = "mps" if torch.backends.mps.is_available() else "cpu"
model = FourierGaussianPairs(centers, coeff, L).to(device)
K = torch.tensor(k[train].reshape(-1, 3), device=device)
Y = torch.tensor(
    np.stack([y[train].real, y[train].imag], -1).reshape(-1, 2), device=device
)
C = torch.tensor(c[train].ravel(), device=device)
vk = torch.tensor(k[val].reshape(-1, 3), device=device)
vy = torch.tensor(
    np.stack([y[val].real, y[val].imag], -1).reshape(-1, 2), device=device
)
vc = torch.tensor(c[val].ravel(), device=device)
optim = torch.optim.Adam(
    [
        {"params": [model.coefficients], "lr": 0.015},
        {"params": [model.centers, model.log_diagonal, model.lower], "lr": 0.002},
    ]
)
center0 = model.centers.detach().clone()
log0 = model.log_diagonal.detach().clone()


def evaluate():
    with torch.no_grad():
        total = 0.0
        den = 0.0
        for start in range(0, len(vk), 1024):
            pred = model(vk[start : start + 1024]) * vc[start : start + 1024, None]
            total += (pred - vy[start : start + 1024]).square().sum().item()
            den += vy[start : start + 1024].square().sum().item()
    return total / den


history = [{"step": 0, "validation_nmse": evaluate()}]
best = history[0]["validation_nmse"]
beststate = None
t = time.time()
for step in range(1, args.steps + 1):
    idx = torch.randint(len(K), (1024,), device=device)
    pred = model(K[idx]) * C[idx, None]
    loss = (
        (pred - Y[idx]).square().mean()
        + 0.0001 * model.coefficients.square().mean()
        + 0.001 * (model.centers - center0).square().mean()
        + 0.001 * (model.log_diagonal - log0).square().mean()
    )
    optim.zero_grad()
    loss.backward()
    optim.step()
    with torch.no_grad():
        model.centers.copy_(
            torch.maximum(torch.minimum(model.centers, center0 + 0.4), center0 - 0.4)
        )
        model.log_diagonal.clamp_(float(np.log(1 / 0.95)), float(np.log(1 / 0.4)))
        model.lower.clamp_(-0.5, 0.5)
    if step % 100 == 0 or step == args.steps:
        score = evaluate()
        history.append(
            {"step": step, "validation_nmse": score, "training_batch_loss": loss.item()}
        )
        print(history[-1], flush=True)
        if score < best:
            best = score
            beststate = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
output = ROOT / "results/adaptive-demo" / args.dataset
output.mkdir(parents=True, exist_ok=True)
if beststate is not None:
    torch.save(beststate, output / "weights.pt")
result = {
    "dataset": args.dataset,
    "device": device,
    "box": 16,
    "spectral_radius": 6,
    "gaussian_pairs": len(centers),
    "training_particles": len(train),
    "validation_particles": len(val),
    "steps": args.steps,
    "seconds": time.time() - t,
    "best_validation_nmse": best,
    "initial_validation_nmse": history[0]["validation_nmse"],
    "history": history,
    "scope": "single-volume, low-frequency nonlinear real-data feasibility; not a three-dataset adaptive benchmark or resolution claim",
}
(output / "metrics.json").write_text(json.dumps(result, indent=2))
print("COMPLETE", result["seconds"], best, flush=True)
