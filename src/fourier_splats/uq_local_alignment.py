"""Batched local alignment against a fixed continuous constant-cell template.

The caller chooses the template independently or labels oracle/misspecified
controls. This is local refinement, never ab initio orientation discovery.
No noise-dependent alignment output carries a confidence guarantee by itself.
"""
import numpy as np
import finufft
from scipy.spatial.transform import Rotation
from .uq_cell_moments import normalized_cell_moments


class CachedCellTemplate:
    def __init__(self, coefficients, box, eps=1e-12, nthreads=1):
        rho = np.asarray(coefficients, float)
        if box < 2 or box % 2 or rho.size != box**3 or not np.isfinite(rho).all():
            raise ValueError('Finite coefficients on an even constant-cell grid required')
        if not 0 < eps < 1 or nthreads < 1:
            raise ValueError('Positive transform accuracy and thread count required')
        coords = (np.arange(box)+.5)/box-.5
        z, y, x = np.meshgrid(coords, coords, coords, indexing='ij')
        self.strengths = np.stack([rho.reshape((box,)*3)*m for m in [np.ones_like(x), x, y, z]]).astype(complex)
        self.box, self.eps, self.nthreads = box, eps, nthreads

    def fourier(self, k, gradient=False):
        k = np.asarray(k, float)
        if k.shape[-1] != 3 or not np.isfinite(k).all():
            raise ValueError('Finite three-dimensional frequencies required')
        b = self.box
        frequencies = [np.ascontiguousarray(2*np.pi*k.reshape(-1, 3)[:, axis]/b) for axis in [2, 1, 0]]
        strengths = self.strengths if gradient else self.strengths[:1]
        transformed = finufft.nufft3d2(*frequencies, strengths, isign=-1, eps=self.eps,
            nthreads=self.nthreads).T.reshape((*k.shape[:-1], len(strengths)))
        transformed *= np.exp(-1j*np.pi*k.sum(axis=-1)/b)[..., None]/b**1.5
        local = normalized_cell_moments(k, b).conj()
        mass = np.prod(local[..., 0], axis=-1)
        value = transformed[..., 0]*mass
        if not gradient:
            return value
        moments = np.empty((*k.shape[:-1], 3), complex)
        for axis in range(3):
            other = np.prod(np.delete(local[..., 0], axis, axis=-1), axis=-1)
            moments[..., axis] = transformed[..., axis+1]*mass+transformed[..., 0]*local[..., axis, 1]*other
        return value, -2j*np.pi*moments

    def signal(self, k, q, ctf, shift, noise, jacobian=False):
        k, q, ctf, shift = (np.asarray(v, float) for v in (k, q, ctf, shift))
        if (ctf.ndim != 2 or k.shape != (*ctf.shape, 3) or q.shape != (*ctf.shape, 2)
                or shift.shape != (len(k), 2) or not np.isfinite(noise) or noise <= 0
                or not all(np.isfinite(v).all() for v in (q, ctf, shift))):
            raise ValueError('Finite consistent particle geometry and positive noise required')
        phase = np.exp(-2j*np.pi*np.einsum('nqa,na->nq', q, shift))*ctf/noise
        if jacobian:
            value, gradient = self.fourier(k, gradient=True)
        else:
            value = self.fourier(k)
        signal = value*phase
        real_signal = np.concatenate([signal.real, signal.imag], axis=1)
        if not jacobian:
            return real_signal
        rotations = np.stack([np.sum(gradient*np.cross(k, np.eye(3)[a]), axis=-1) for a in range(3)], axis=-1)
        translations = -2j*np.pi*q*value[..., None]
        j = np.concatenate([rotations, translations], axis=-1)*phase[..., None]
        return real_signal, np.concatenate([j.real, j.imag], axis=1)


def refine_local_poses(template, k, q, ctf, observations, noise, field_A, *,
                        initial_rotation, initial_shift_A, iterations=10,
                        damping=.01, step_rotation_degrees=2., step_shift_A=1.,
                        search_rotation_degrees=20., search_shift_A=10.):
    """Fixed-iteration damped Gauss--Newton with per-particle backtracking.

    Refinement is constrained relative to the supplied initialization, not to
    truth. Rotation updates compose on the right. Every line-search candidate
    is evaluated, including the unchanged iterate. No outcome-based retries.
    """
    k, q, ctf, y = (np.asarray(v, float) for v in (k, q, ctf, observations))
    n, nq = ctf.shape
    r0 = np.asarray(initial_rotation, float)
    t0 = np.asarray(initial_shift_A, float)
    if (k.shape != (n, nq, 3) or q.shape != (n, nq, 2) or y.shape != (n, 2*nq)
            or r0.shape != (n, 3, 3) or t0.shape != (n, 2)
            or not all(np.isfinite(v).all() for v in (k, q, ctf, y, r0, t0))
            or iterations < 1 or not all(np.isfinite(v) and v > 0 for v in
                (noise, field_A, damping, step_rotation_degrees, step_shift_A,
                 search_rotation_degrees, search_shift_A))):
        raise ValueError('Finite compatible local-alignment inputs required')
    if not np.allclose(r0.swapaxes(-1, -2)@r0, np.eye(3), atol=1e-10) or not np.allclose(np.linalg.det(r0), 1., atol=1e-10):
        raise ValueError('Initial rotations must lie in SO(3)')
    rotations, shifts = r0.copy(), t0.copy()
    scale = np.array([np.pi/180]*3+[1/field_A]*2)
    history = []
    def signal(r, t, jac=False):
        return template.signal(k@r, q, ctf, t/field_A, noise, jacobian=jac)
    current = signal(rotations, shifts)
    objective = np.sum((y-current)**2, axis=1)
    initial_objective = objective.copy()
    for it in range(iterations):
        current, j = signal(rotations, shifts, True)
        j = j*scale
        information = j.swapaxes(-1, -2)@j
        gradient = np.einsum('nmp,nm->np', j, y-current)
        regularizer = damping*np.maximum(np.diagonal(information, axis1=-2, axis2=-1), 1e-8)
        step = np.linalg.solve(information+regularizer[:, :, None]*np.eye(5), gradient[..., None])[..., 0]
        for block, maximum in [(slice(0, 3), step_rotation_degrees), (slice(3, 5), step_shift_A)]:
            step[:, block] *= np.minimum(1., maximum/np.maximum(np.linalg.norm(step[:, block], axis=1), 1e-300))[:, None]
        best_r, best_t, best_objective = rotations.copy(), shifts.copy(), objective.copy()
        accepted = np.zeros(n)
        for factor in [1., .5, .25]:
            candidate_r = rotations@Rotation.from_rotvec(np.deg2rad(factor*step[:, :3])).as_matrix()
            relative = Rotation.from_matrix(r0.swapaxes(-1, -2)@candidate_r).as_rotvec()
            relative *= np.minimum(1., np.deg2rad(search_rotation_degrees)/np.maximum(np.linalg.norm(relative, axis=1), 1e-300))[:, None]
            candidate_r = r0@Rotation.from_rotvec(relative).as_matrix()
            delta_t = shifts+factor*step[:, 3:]-t0
            delta_t *= np.minimum(1., search_shift_A/np.maximum(np.linalg.norm(delta_t, axis=1), 1e-300))[:, None]
            candidate_t = t0+delta_t
            pred = signal(candidate_r, candidate_t)
            candidate_objective = np.sum((y-pred)**2, axis=1)
            improve = candidate_objective < best_objective
            best_r[improve], best_t[improve] = candidate_r[improve], candidate_t[improve]
            best_objective[improve] = candidate_objective[improve]
            accepted[improve] = factor
        history.append(dict(iteration=it+1, objective_sum=float(best_objective.sum()),
            improved_particles=int(np.count_nonzero(accepted)),
            accepted_factors=accepted.tolist()))
        rotations, shifts, objective = best_r, best_t, best_objective
    relative = Rotation.from_matrix(r0.swapaxes(-1, -2)@rotations).magnitude()*180/np.pi
    shift_move = np.linalg.norm(shifts-t0, axis=1)
    return dict(rotations=rotations, shifts_A=shifts,
        initial_objective=initial_objective, final_objective=objective, history=history,
        rotation_boundary_particles=int(np.sum(relative >= search_rotation_degrees*(1-1e-6))),
        shift_boundary_particles=int(np.sum(shift_move >= search_shift_A*(1-1e-6))))
