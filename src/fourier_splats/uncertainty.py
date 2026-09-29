"""Bias-aware linear-functional inference with particle-specific nuisance blocks.

The input experiment is whitened and real-valued. Bounds are assumptions supplied
by the caller, not inferred by this module. All design choices, including the
pilot and noise whitening, must be independent of the inference noise for the
stated fixed-design Gaussian coverage result. General support-function inference
and bias-aware normal critical values are classical optimal-recovery tools.
"""

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import issparse
from scipy.sparse.linalg import LinearOperator, cg
from scipy.stats import norm


def bias_aware_half_width(sd, bias, alpha=0.05):
    """Exact normal half-width uniformly over absolute bias <= bias.

    For |b| <= bias and Z~N(0,sd^2), P(|Z+b| <= q) >= 1-alpha.
    This is a folded-normal critical value, not a posterior credible interval.
    """
    sd, bias = float(sd), float(bias)
    if sd < 0 or bias < 0 or not 0 < alpha < 1:
        raise ValueError("Nonnegative sd/bias and alpha in (0,1) required")
    if sd == 0:
        return bias
    t = bias / sd
    # For a very large bias, the remote Gaussian tail is numerically negligible.
    if t > 40:
        return bias + norm.ppf(1 - alpha) * sd
    def coverage(q):
        return norm.cdf(q - t) - norm.cdf(-q - t) - (1 - alpha)
    return sd * brentq(coverage, 0, t + norm.ppf(1 - alpha / 2) + 1)


@dataclass
class Certificate:
    weights: np.ndarray
    noise_sd: float
    reconstruction_bias: float
    nuisance_bias: float
    remainder_bias: float
    alpha: float
    objective: float
    dual_lower_bound: float
    iterations: int
    converged: bool
    history: list
    nuisance_group_biases: list | None = None

    @property
    def bias(self):
        return self.reconstruction_bias + self.nuisance_bias + self.remainder_bias

    @property
    def half_width(self):
        return bias_aware_half_width(self.noise_sd, self.bias, self.alpha)

    def estimate(self, residual, pilot_target=0.0):
        return float(pilot_target + self.weights @ np.asarray(residual).ravel())

    def interval(self, residual, pilot_target=0.0):
        center = self.estimate(residual, pilot_target)
        return center - self.half_width, center + self.half_width


def _broadcast_bounds(x, n, name):
    out = np.broadcast_to(np.asarray(x, dtype=float), (n,)).copy()
    if not np.isfinite(out).all() or (out < 0).any():
        raise ValueError(name + " must contain finite nonnegative bounds")
    return out


def compress_nuisance_group(design,relative_jitter=1e-12):
    """Replace a wide group design by its pixel-space Gram factor.

    ||D' w|| equals sqrt(w'DD'w). Cholesky of DD'+jitter*I slightly enlarges
    this support function, preserving a conservative bound. This exact Gram
    reduction (apart from the positive numerical padding) avoids storing the
    number of density coefficients times the number of nuisance directions.
    """
    design=np.asarray(design,dtype=float)
    if design.ndim!=3 or relative_jitter<=0:
        raise ValueError('A three-dimensional design and positive padding are required')
    if design.shape[2]<=design.shape[1]:
        return design.copy()
    gram=design@design.transpose(0,2,1)
    scale=np.maximum(np.max(np.diagonal(gram,axis1=1,axis2=2),axis=1),1e-300)
    return np.linalg.cholesky(gram+relative_jitter*scale[:,None,None]*np.eye(design.shape[1])[None])


def optimize_certificate(
    a, j, ell, density_radius, nuisance_radius=0.0, remainder_radius=0.0,
    alpha=0.05, maxiter=100, rtol=1e-4, smoothing=1e-7,
    cg_rtol=1e-10, cg_maxiter=1000,
    extra_nuisance_groups=(),
):
    """Minimize a valid sum-of-norms half-width by quadratic majorization.

    a: (n*m,p) design, dense or sparse; j: (n,m,q) nuisance Jacobians.
    ell: (p,) functional; density_radius: bound on ||c-c0||_2.
    nuisance_radius and remainder_radius: scalar or (n,) Euclidean bounds.
    For a general ellipsoid, transform a and ell before calling this routine.

    The optimized objective is z*sd+bias. The returned interval uses the tighter
    exact bias-aware Gaussian critical value. The primal-dual gap certifies
    optimization accuracy only for the sum-of-norms objective. Coverage does not
    require optimization convergence, only valid bounds and fixed weights.
    extra_nuisance_groups is a sequence of (design, radius) pairs. Each design
    has shape (n,m,q_g), and its per-particle radius is scalar or (n,).
    These independent norm groups can relax structured nonlinear interactions.
    """
    j = np.asarray(j, dtype=np.float64)
    ell = np.asarray(ell, dtype=np.float64)
    if j.ndim != 3 or a.shape[0] != j.shape[0] * j.shape[1] or a.shape[1] != len(ell):
        raise ValueError("Incompatible design, functional, and particle blocks")
    n, m, q = j.shape
    B = float(density_radius)
    if not np.isfinite(B) or B < 0 or not 0 < alpha < 1:
        raise ValueError("A finite nonnegative density bound and valid alpha are required")
    eta = _broadcast_bounds(nuisance_radius, n, "nuisance_radius")
    designs=[j]; radii=[eta]; slices=[slice(0,q)]
    for design,radius in extra_nuisance_groups:
        design=np.asarray(design,dtype=np.float64)
        if design.ndim!=3 or design.shape[:2]!=(n,m):
            raise ValueError('Each extra nuisance design must have shape (n,m,q_g)')
        designs.append(design)
        radii.append(_broadcast_bounds(radius,n,'extra nuisance radius'))
        slices.append(slice(q,q+design.shape[2]));q+=design.shape[2]
    j=np.concatenate(designs,axis=-1)
    eta=np.stack(radii,axis=1)
    gamma = _broadcast_bounds(remainder_radius, n, "remainder_radius")
    z = norm.ppf(1 - alpha / 2)
    at = a.T
    jtj = np.einsum("nmq,nmr->nqr", j, j) if q<=m else None
    diag_a = np.asarray(a.power(2).sum(axis=0)).ravel() if issparse(a) else np.sum(np.asarray(a)**2, axis=0)
    base = B * np.linalg.norm(ell)
    epsilon = max(base, 1e-12) * smoothing
    w = np.zeros(a.shape[0], dtype=np.float64)
    history = []

    def pieces(v):
        vb = v.reshape(n, m)
        residual = ell - at @ v
        jt = np.einsum("nmq,nm->nq", j, vb)
        sd = np.linalg.norm(v)
        b = B * np.linalg.norm(residual)
        nb = eta * np.stack([np.linalg.norm(jt[:,s],axis=1) for s in slices],axis=1)
        rb = gamma * np.linalg.norm(vb, axis=1)
        return residual, jt, sd, b, nb, rb

    def dual_bound(v, residual, jt):
        # A feasible dual point gives a numerical lower bound on the optimum.
        rd = np.linalg.norm(residual)
        # Use the smooth objective's feasible subgradient. Normalizing a nearly
        # zero residual to norm B picks an arbitrary boundary subgradient and
        # needlessly destroys the dual certificate at exact-fit optima.
        dual_v = B*B * residual / max(np.sqrt((B*rd)**2 + epsilon**2), 1e-300)
        dual_t=np.empty_like(jt)
        for g,s in enumerate(slices):
            jnorm=np.linalg.norm(jt[:,s],axis=1)
            dual_t[:,s]=jt[:,s]*(eta[:,g]**2/np.maximum(np.sqrt((eta[:,g]*jnorm)**2+epsilon**2),1e-300))[:,None]
        vb = v.reshape(n, m)
        vnorm = np.linalg.norm(vb, axis=1)
        dual_r = vb * (gamma*gamma / np.maximum(np.sqrt((gamma*vnorm)**2 + epsilon**2), 1e-300))[:, None]
        leftover = (a @ dual_v).reshape(n, m) - np.einsum("nmq,nq->nm", j, dual_t) - dual_r
        factor = min(1.0, z / max(np.linalg.norm(leftover), 1e-300))
        return max(0.0, float(factor * (ell @ dual_v)))

    best_w = w.copy()
    best_objective = base
    best_lower = 0.0
    converged = base == 0
    coef = np.zeros(a.shape[1])
    for it in range(maxiter if base else 0):
        residual, jt, sd, b, nb, rb = pieces(w)
        objective = z * sd + b + nb.sum() + rb.sum()
        lower = dual_bound(w, residual, jt)
        if objective <= best_objective:
            best_objective, best_w = float(objective), w.copy()
        best_lower = max(best_lower, lower)
        history.append({"iteration": it, "objective": float(objective), "lower_bound": lower, "best_gap": float(max(0, best_objective - best_lower))})
        if best_objective - best_lower <= rtol * max(best_objective, 1e-15):
            converged = True
            break
        # Smooth norm majorization. Bounds equal to zero contribute no penalty.
        t0 = np.sqrt((z * sd)**2 + epsilon**2)
        tb = np.sqrt(b*b + epsilon**2)
        tn = np.sqrt(nb*nb + epsilon**2)
        tr = np.sqrt(rb*rb + epsilon**2)
        diagonal = z*z / t0 + gamma*gamma / tr
        group_weight = eta*eta / tn
        nuisance_weight=np.concatenate([np.broadcast_to(group_weight[:,g,None],(n,s.stop-s.start)) for g,s in enumerate(slices)],axis=1)
        beta = B*B / tb
        root_weight=np.sqrt(nuisance_weight)
        if q<=m:
            inv_small=np.linalg.inv(np.eye(q)[None]+jtj*root_weight[:,:,None]*root_weight[:,None,:]/diagonal[:,None,None])
        else:
            # For many derivative directions the pixel block is smaller than
            # the Woodbury matrix. Solve that positive definite block directly.
            weighted_j=j*root_weight[:,None,:]
            d=diagonal[:,None,None]*np.eye(m)[None]+weighted_j@weighted_j.transpose(0,2,1)
            inv_pixel=np.linalg.inv(d)

        def dinv(x):
            xb = np.asarray(x).reshape(n, m)
            if q>m:
                return np.einsum('nmk,nk->nm',inv_pixel,xb).ravel()
            small = np.einsum("nmq,nm->nq", j, xb)*root_weight
            small = np.einsum("nqr,nr->nq", inv_small, small)
            correction = np.einsum("nmq,nq->nm", j, small*root_weight)
            result = xb / diagonal[:, None] - correction/diagonal[:,None]**2
            return result.ravel()

        ridge = 1.0 / beta
        op = LinearOperator((a.shape[1],)*2, matvec=lambda x: at @ dinv(a @ x) + ridge*x, dtype=np.float64)
        # Conservative diagonal approximation; ignored low-rank correction affects
        # speed, not the matrix being solved or the returned certificate.
        pre_diag = diag_a / np.median(diagonal) + ridge
        pre = LinearOperator(op.shape, matvec=lambda x: x / np.maximum(pre_diag, 1e-300), dtype=np.float64)
        coef, info = cg(op, ell, x0=coef, rtol=cg_rtol, atol=0, maxiter=cg_maxiter, M=pre)
        normal_residual = np.linalg.norm(op @ coef - ell) / max(np.linalg.norm(ell), 1e-300)
        history[-1].update(cg_info=int(info), cg_relative_residual=float(normal_residual))
        w = dinv(a @ coef)
        # The coefficient-space solution is a stable dual candidate even when
        # ell-A'w is close to zero. Forming that tiny residual and normalizing it
        # loses the interior subgradient of the nonsmooth bias norm.
        wb = w.reshape(n, m)
        dual_t = nuisance_weight * np.einsum("nmq,nm->nq", j, wb)
        dual_r = (gamma*gamma/tr)[:, None] * wb
        dual_u = (a @ coef).reshape(n,m) - np.einsum("nmq,nq->nm", j, dual_t) - dual_r
        nt = np.stack([np.linalg.norm(dual_t[:,s],axis=1) for s in slices],axis=1)
        nr = np.linalg.norm(dual_r, axis=1)
        ratio_t = np.divide(eta, nt, out=np.full_like(eta, np.inf), where=nt>0)
        ratio_r = np.divide(gamma, nr, out=np.full(n, np.inf), where=nr>0)
        dual_scale = min(1.0, B/max(np.linalg.norm(coef),1e-300), z/max(np.linalg.norm(dual_u),1e-300), float(ratio_t.min()), float(ratio_r.min()))
        best_lower = max(best_lower, float(max(0, dual_scale*(ell @ coef))))

    residual, jt, sd, b, nb, rb = pieces(w)
    objective = z*sd + b + nb.sum() + rb.sum()
    if objective < best_objective:
        best_objective, best_w = float(objective), w.copy()
    best_lower = max(best_lower, dual_bound(w, residual, jt))
    residual, jt, sd, b, nb, rb = pieces(best_w)
    return Certificate(best_w, float(sd), float(b), float(nb.sum()), float(rb.sum()), float(alpha), best_objective, min(best_lower, best_objective), len(history), bool(converged or best_objective - best_lower <= rtol*max(best_objective, 1e-15)), history,nb.sum(axis=0).tolist())
