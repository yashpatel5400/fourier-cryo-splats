import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_continuous_gaussian import cell_pose_jacobian
from fourier_splats.uq_local_alignment import CachedCellTemplate, refine_local_poses


def test_cached_template_matches_independent_continuous_implementations():
    rng = np.random.default_rng(301)
    k = rng.normal(size=(3, 10, 3)); q = rng.normal(size=(3, 10, 2))
    ctf = rng.normal(size=(3, 10)); rho = rng.normal(size=8**3)
    model = CachedCellTemplate(rho, 8)
    signal, jacobian = model.signal(k, q, ctf, np.zeros((3, 2)), .7, jacobian=True)
    np.testing.assert_allclose(signal.ravel(), cell_forward(k, ctf, rho, 8, .7), atol=1e-10)
    np.testing.assert_allclose(jacobian, cell_pose_jacobian(k, q, ctf, rho, 8, .7), atol=1e-9)
    shift = rng.normal(size=(3, 2))*.03
    _, j = model.signal(k, q, ctf, shift, .7, jacobian=True)
    for axis in range(5):
        h = 1e-6; u = np.zeros((3, 3)); t = np.zeros((3, 2))
        if axis < 3:
            u[:, axis] = h
        else:
            t[:, axis-3] = h
        plus = model.signal(k@Rotation.from_rotvec(u).as_matrix(), q, ctf, shift+t, .7)
        minus = model.signal(k@Rotation.from_rotvec(-u).as_matrix(), q, ctf, shift-t, .7)
        np.testing.assert_allclose(j[..., axis], (plus-minus)/(2*h), rtol=1e-6, atol=3e-8)


def test_noise_free_local_refinement_and_search_bounds():
    rng = np.random.default_rng(302)
    rho = rng.normal(size=12**3)
    k = rng.normal(size=(4, 30, 3)); q = rng.normal(size=(4, 30, 2)); ctf = np.ones((4, 30))
    model = CachedCellTemplate(rho, 12)
    y = model.signal(k, q, ctf, np.zeros((4, 2)), .1)
    r0 = Rotation.from_rotvec(rng.normal(size=(4, 3))*.02).as_matrix()
    t0 = rng.normal(size=(4, 2))*.4
    fit = refine_local_poses(model, k, q, ctf, y, .1, 200.,
        initial_rotation=r0, initial_shift_A=t0, iterations=8)
    assert np.max(Rotation.from_matrix(fit['rotations']).magnitude()) < 1e-7
    assert np.max(abs(fit['shifts_A'])) < 1e-5
    sums = [fit['initial_objective'].sum()]+[r['objective_sum'] for r in fit['history']]
    assert np.all(np.diff(sums) <= 1e-12)
    bounded = refine_local_poses(model, k, q, ctf, y, .1, 200., initial_rotation=r0,
        initial_shift_A=t0, iterations=3, search_rotation_degrees=.01, search_shift_A=.01)
    assert np.max(Rotation.from_matrix(r0.swapaxes(-1, -2)@bounded['rotations']).magnitude()) <= np.deg2rad(.01)*(1+1e-8)
    assert np.max(np.linalg.norm(bounded['shifts_A']-t0, axis=1)) <= .01*(1+1e-8)
