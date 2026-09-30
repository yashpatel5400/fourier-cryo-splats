import numpy as np
from fourier_splats.uq_support import (rescale_supported_problem,
    restrict_cell_density, outside_tail_bias_upper)
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_pose import pose_cell_forward, cube_quadrature, perturbed_geometry


def test_support_isometry_preserves_cell_signal_target_and_nonlinear_pose():
    rng = np.random.default_rng(921)
    k = rng.normal(size=(3, 4, 3)); q = rng.normal(size=(3, 4, 2))
    ctf = rng.normal(size=(3, 4)); side = .5; box = 8
    a = np.zeros((box,)*3); a[2:6, 2:6, 2:6] = rng.normal(size=(4,)*3)
    small, count, tail = restrict_cell_density(a, box, side)
    assert tail == 0
    np.testing.assert_allclose(np.linalg.norm(small), np.linalg.norm(a))
    centers = [[0, 0, .03], [0, 0, -.03]]; signs = [1., -1.]
    p = rescale_supported_problem(k, q, ctf, centers, signs, .08, .002, side)
    np.testing.assert_allclose(cell_forward(k, ctf, a, box, .2),
        cell_forward(p['k'], p['ctf'], small, count, .2), rtol=2e-9, atol=2e-10)
    np.testing.assert_allclose(cell_target_coefficients(box, centers, signs, .08)@a.ravel(),
        cell_target_coefficients(count, p['centers'], p['signs'], p['width'])@small,
        rtol=2e-12)
    u = rng.normal(size=(3, 5)); u /= np.linalg.norm(u, axis=1)[:, None]
    np.testing.assert_allclose(pose_cell_forward(k, q, ctf, a, box, .2, u, .04, .002),
        pose_cell_forward(p['k'], p['q'], p['ctf'], small, count, .2, u, .04, p['shift']),
        rtol=2e-9, atol=2e-10)


def test_tail_allowance_dominates_independent_nonlinear_fourier_fields():
    rng = np.random.default_rng(922)
    n, nq = 4, 5
    k = rng.normal(size=(n, nq, 3))*3; q = rng.normal(size=(n, nq, 2))*3
    ctf = rng.normal(size=(n, nq)); w = rng.normal(size=n*2*nq)
    xyz, quad = cube_quadrature(28); noise = .4; angle = .15; shift = .01
    bound = outside_tail_bias_upper(k, q, ctf, w, noise, angle, shift, 0., 1.)
    ww = w.reshape(n, 2*nq); c = (ww[:, :nq]+1j*ww[:, nq:])*ctf/noise
    nominal = np.real(np.sum(c[..., None]*np.exp(2j*np.pi*np.einsum('nqj,pj->nqp', k, xyz)), axis=(0, 1)))
    for _ in range(6):
        u = rng.normal(size=(n, 5)); u /= np.linalg.norm(u, axis=1)[:, None]
        ku, shifts = perturbed_geometry(k, q, u, angle, shift)
        actual = np.real(np.sum(c[..., None]*np.exp(2j*np.pi*(
            np.einsum('nqj,pj->nqp', ku, xyz)+shifts[..., None])), axis=(0, 1)))
        assert np.sqrt(quad@((actual-nominal)**2)) <= bound['pose_field_upper']
