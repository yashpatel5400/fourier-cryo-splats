import numpy as np
import pytest
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.linalg import block_diag
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_cubic_design import DifferentiableCubicPoseFieldOperator
from fourier_splats.uq_cubic_pose import gaussian_moments_cubic
from fourier_splats.uq_sobolev_penalty import SobolevRemainderPenalty
from fourier_splats.uq_cubic_optimization import CubicDesignObjective, optimize_cubic_weights


def fixture(smoothing=0., noise=.1):
    q = np.array([[[.31, .06], [-.04, .37]]]); k = np.pad(q, ((0, 0), (0, 0), (0, 1)))
    ctf = np.array([[1., -.8]]); w = np.array([.35, -.2, .15, -.25])
    gram = QuadratureObservationGram(k, ctf, noise, order=12, preconditioner_rank=0)
    op = DifferentiableCubicPoseFieldOperator(k, q, ctf, w, noise, .17, .08, order=4, backend='direct')
    scaling = op.establish_quadrature_design_scaling(4)
    remainder = SobolevRemainderPenalty(k, q, ctf/noise, .17, .08)
    moments = .1*gaussian_moments_cubic(k, [[.03, -.04, .02]], [1.], .2)
    obj = CubicDesignObjective(gram, op, remainder, moments, [[.08, -.06, 0.]], [1.], .24,
        2., 1., np.sqrt(scaling['sum_group_scales']), norm.isf((.05/12-1e-6/12)/2), smoothing=smoothing)
    return obj, w


def full_modes(obj, w):
    obj.op.set_weights(w)
    matrix = np.column_stack([obj.op.matvec(e) for e in np.eye(obj.op.shape[1])])
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    return u, s


@pytest.mark.parametrize('smoothing', [0., .002])
def test_full_objective_gradient_and_norm_supports(smoothing):
    obj, w = fixture(smoothing); rng = np.random.default_rng(650101)
    def evaluate(v):
        u, _ = full_modes(obj, v)
        return obj.evaluate(v, u[:, :3])
    value, gradient, _, support = evaluate(w)
    direction = rng.normal(size=len(w)); step = 1e-5
    fd = (evaluate(w+step*direction)[0]-evaluate(w-step*direction)[0])/(2*step)
    np.testing.assert_allclose(fd, gradient@direction, rtol=2e-7, atol=1e-9)
    for _ in range(5):
        other = rng.normal(size=len(w)); v, _, row, _ = evaluate(other)
        assert v >= value+gradient@(other-w)-1e-9
        # Support must apply to the UNSMOOTHED norm sum even if design smooths.
        old = obj.smoothing; obj.smoothing = 0.
        _, _, exact, _ = evaluate(other); obj.smoothing = old
        norm_sum = exact['pose_polynomial_term']+exact['pilot_polynomial_term']+exact['quartic_term']
        assert support@other <= norm_sum+1e-10


def matrix_sqrt(g):
    e, v = np.linalg.eigh(.5*(g+g.T))
    assert e.min() >= -1e-10*max(1., e.max())
    return np.sqrt(np.maximum(e, 0.))[:, None]*v.T


def independent_conic(obj):
    import cvxpy as cp
    dimension = len(obj.a); basis = np.eye(dimension)
    # Analytic continuous cube Gram, independently assembled in real coordinates.
    k = obj.op.k[0]; t = obj.op.transfer[0]
    minus = np.prod(np.sinc(k[:, None]-k[None, :]), axis=-1)
    plus = np.prod(np.sinc(k[:, None]+k[None, :]), axis=-1)
    gram = block_diag(.5*t[:, None]*t[None, :]*(minus+plus),
                      .5*t[:, None]*t[None, :]*(minus-plus))
    for i in range(dimension):
        np.testing.assert_allclose(obj.gram.matvec(basis[i]), gram[:, i], atol=3e-12)
    augmented = np.block([[gram+obj.eta*np.eye(dimension), -obj.a[:, None]],
                          [-obj.a[None], np.array([[obj.target_norm2]])]])
    density_root = matrix_sqrt(augmented)
    matrices, pilot_columns = [], []
    for e in basis:
        obj.op.set_weights(e)
        matrices.append(np.column_stack([obj.op.matvec(v) for v in np.eye(obj.op.shape[1])]))
        pilot_columns.append(obj.op.pair_moments(obj.moments))
    pilot = np.column_stack(pilot_columns)
    # Compress only common numerical null spaces; retain >1e-13 relative modes.
    left, sv, _ = np.linalg.svd(np.concatenate(matrices, axis=1), full_matrices=False)
    right, sr, _ = np.linalg.svd(np.concatenate([m.T for m in matrices], axis=1), full_matrices=False)
    left = left[:, sv > sv[0]*1e-13]; right = right[:, sr > sr[0]*1e-13]
    compressed = [left.T@m@right for m in matrices]
    assert max(np.linalg.norm(m-left@c@right.T) for m, c in zip(matrices, compressed)) < 1e-12
    w = cp.Variable(dimension)
    amplitude = cp.hstack([cp.norm(cp.hstack([w[j], w[obj.nq+j]])) for j in range(obj.nq)])
    field = sum(w[j]*compressed[j] for j in range(dimension))
    masses = obj.coefficients[0]@amplitude
    pose = cp.norm(cp.hstack([cp.norm(field, 2), np.sqrt(obj.kernel_error)*masses]))
    remainder = 0.
    for j in range(obj.remainder.order):
        root = matrix_sqrt(block_diag(*obj.remainder.blocks[0, j]))
        remainder += obj.remainder.multipliers[0, j]*cp.norm(root@w)
    remainder += obj.remainder.residual[0]@amplitude
    cost = (obj.z*cp.norm(w)+obj.B*cp.norm(density_root@cp.hstack([w, 1.]))+
            obj.B*obj.L*pose+obj.L*cp.norm(pilot@w)+(obj.B+obj.P)*remainder)
    problem = cp.Problem(cp.Minimize(cost))
    optimum = problem.solve(solver='CLARABEL', tol_gap_abs=1e-9, tol_gap_rel=1e-9, tol_feas=1e-9)
    assert problem.status == 'optimal'
    return optimum, w.value


def test_independent_conic_optimum_and_continuous_dual_bracket():
    obj, initial = fixture(); optimum, conic_weights = independent_conic(obj)
    def evaluate(w):
        modes, _ = full_modes(obj, w)
        value, gradient, _, _ = obj.evaluate(w, modes[:, :1])
        return value, gradient
    fitted = minimize(evaluate, initial, method='L-BFGS-B', jac=True,
        options={'maxiter': 200, 'ftol': 1e-13, 'gtol': 1e-9})
    modes, singular = full_modes(obj, fitted.x)
    value, _, _, support = obj.evaluate(fitted.x, modes[:, :1])
    audit = obj.audit(fitted.x, support, singular[0]**2)
    assert audit['dual_lower_bound'] <= optimum+2e-7
    assert audit['sum_objective_upper'] >= optimum-2e-7
    assert abs(value-optimum) < 2e-6
    assert audit['relative_sum_gap'] < 2e-5
    print('CONIC_VALIDATION', {'conic_optimum': optimum, 'design_value': value, 'dual_lower': audit['dual_lower_bound'], 'audited_upper': audit['sum_objective_upper'], 'relative_gap': audit['relative_sum_gap'], 'optimizer_success': bool(fitted.success)})
    np.testing.assert_allclose(evaluate(conic_weights)[0], optimum, atol=2e-7)


def test_partial_ritz_modes_still_use_independent_final_upper(monkeypatch):
    import fourier_splats.uq_cubic_optimization as implementation
    from scipy.sparse.linalg import ArpackNoConvergence
    obj, initial = fixture(.002)
    def partial(operator, **kwargs):
        v = np.ones(operator.shape[0]); v /= np.linalg.norm(v)
        raise ArpackNoConvergence('deliberate partial modes', np.array([v@(operator@v)]), v[:, None])
    monkeypatch.setattr(implementation, 'eigsh', partial)
    result = optimize_cubic_weights(obj, initial, maxiter=1, max_evaluations=3,
        modes_count=2, ritz_subspace=5, power_iterations=8, certificate_seed=650103, precondition=False)
    _, singular = full_modes(obj, result['weights'])
    assert result['partial_mode_evaluations'] > 0
    assert result['spectral_upper_certificate']['seed'] == 650103
    assert result['spectral_upper_certificate']['eigenvalue_upper'] >= singular[0]**2*(1-1e-10)
    assert result['dual_lower_bound'] <= result['sum_objective_upper']
    assert len(result['optimization_history']) >= result['selected_evaluation']


def test_original_zero_weight_boundary_has_analytic_global_optimum():
    obj, _ = fixture(noise=.7)
    # The first conic attempt reported optimal_inaccurate in this boundary case.
    # Convex density support plus the noise norm instead proves w=0 optimal.
    assert obj.B*np.linalg.norm(obj.a)/np.sqrt(obj.target_norm2) < obj.z
    zero = np.zeros(len(obj.a)); modes, _ = full_modes(obj, zero)
    value, _, _, _ = obj.evaluate(zero, modes[:, :1])
    np.testing.assert_allclose(value, obj.B*np.sqrt(obj.target_norm2), atol=1e-14)
    rng = np.random.default_rng(650107)
    for _ in range(8):
        other = rng.normal(size=len(zero)); modes, _ = full_modes(obj, other)
        assert obj.evaluate(other, modes[:, :1])[0] >= value
