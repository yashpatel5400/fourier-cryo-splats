import numpy as np
from numpy.testing import assert_allclose
from fourier_splats.uq_continuous import fourier_field_norm2, gaussian_target_integrals
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_end_to_end import dephase_observations
from fourier_splats.uq_refitting_diagnostics import dephase_adjoint_weights, realized_class_envelopes


def test_dephase_adjoint_sign_and_noise_norm():
    rng = np.random.default_rng(891)
    q = rng.normal(size=(5, 9, 2)); t = rng.normal(size=(5, 2))
    y = rng.normal(size=(5, 18)); w = rng.normal(size=90)
    v = dephase_adjoint_weights(w, q, t, 13.)
    assert_allclose(w@dephase_observations(y, q, t, 13.), v@y.ravel(), atol=2e-14)
    assert_allclose(v@v, w@w, rtol=3e-15)


def test_realized_envelope_against_independent_dense_sinc():
    rng = np.random.default_rng(919)
    k = rng.normal(size=(2, 7, 3)); khat = k+rng.normal(size=k.shape)*.08
    ctf = rng.uniform(-1, 1, (2, 7)); w = rng.normal(size=28)
    q = rng.normal(size=(2, 7, 2)); shifts = rng.normal(size=(2, 2))
    v = dephase_adjoint_weights(w, q, shifts, 15.)
    gt = QuadratureObservationGram(k, ctf, .7, order=30, preconditioner_rank=0)
    gf = QuadratureObservationGram(khat, ctf, .7, order=30, preconditioner_rank=0)
    centers, signs, width = [[.02, -.01, .06]], [1.], .14
    out = realized_class_envelopes(gt, gf, w, v, centers, signs, width, 2., -.17)
    c = gt.coefficients(v); ch = gf.coefficients(w)
    ff, ell2 = gaussian_target_integrals(k, centers, signs, width)
    norm2, _ = fourier_field_norm2(k, c)
    expected_res2 = ell2-2*np.real(np.sum(c*ff.conj()))+norm2
    union_k = np.concatenate([k.reshape(-1, 3), khat.reshape(-1, 3)])
    union_c = np.concatenate([c.ravel(), -ch.ravel()])
    diff2, _ = fourier_field_norm2(union_k, union_c)
    assert_allclose(out['residual_norm2_unpadded'], expected_res2, atol=2e-9)
    assert_allclose(out['pose_norm2_unpadded'], diff2, atol=2e-9)
    assert_allclose(out['realized_total_bias_envelope'], .17+2*np.sqrt(expected_res2), rtol=1e-6)


def test_identical_pose_envelope_and_sharp_ball_extremum():
    rng = np.random.default_rng(19)
    k = rng.normal(size=(1, 5, 3)); ctf = np.ones((1, 5)); w = rng.normal(size=10)
    g = QuadratureObservationGram(k, ctf, 1., order=30, preconditioner_rank=0)
    d = realized_class_envelopes(g, g, w, w, [[0., 0., 0.]], [1.], .1, 2., 0.)
    assert d['pose_norm2_unpadded'] == 0.
    # The class supremum is attained along the residual, with the sign of delta.
    residual = np.array([1.2, -.9]); B = 2.; delta = -.17
    h = -np.sign(delta)*B*residual/np.linalg.norm(residual)
    assert_allclose(abs(delta-residual@h), abs(delta)+B*np.linalg.norm(residual))
