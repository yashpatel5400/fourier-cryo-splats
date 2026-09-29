import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.physics import plane_restriction, hermitian_pairs, ctf
from fourier_splats.basis import design, evaluate_grid
from fourier_splats.fsc import fsc, resolution


def test_analytic_plane_restriction():
    rng = np.random.default_rng(4)
    L = rng.normal(size=(3, 3))
    P = L @ L.T + np.eye(3)
    E = Rotation.random(random_state=rng).as_matrix()[:, :2]
    mu = rng.normal(size=3)
    q = rng.normal(size=(80, 2))
    m, A, a = plane_restriction(mu, P, E)
    d = q @ E.T - mu
    direct = np.exp(-0.5 * np.einsum("ni,ij,nj->n", d, P, d))
    dd = q - m
    np.testing.assert_allclose(
        direct,
        a * np.exp(-0.5 * np.einsum("ni,ij,nj->n", dd, A, dd)),
        rtol=1e-12,
        atol=1e-14,
    )


def test_conjugate_symmetry():
    rng = np.random.default_rng(2)
    k = rng.normal(size=(50, 3))
    mu = rng.normal(size=(5, 3))
    p = np.tile(np.eye(3), (5, 1, 1))
    c = rng.normal(size=5) + 1j * rng.normal(size=5)
    np.testing.assert_allclose(
        hermitian_pairs(-k, mu, p, c), np.conj(hermitian_pairs(k, mu, p, c)), atol=1e-14
    )


def test_sparse_pairing_and_adjoint():
    rng = np.random.default_rng(3)
    k = rng.uniform(-5, 5, (100, 3))
    for kind in ["gaussian", "voxel"]:
        a, b = design(k, 16, kind)
        am, bm = design(-k, 16, kind)
        np.testing.assert_allclose(a.toarray(), am.toarray(), atol=2e-6)
        np.testing.assert_allclose(b.toarray(), -bm.toarray(), atol=2e-6)
        x = rng.normal(size=a.shape[1])
        y = rng.normal(size=len(k))
        np.testing.assert_allclose(
            np.dot(a @ x, y), np.dot(x, a.T @ y), rtol=1e-12, atol=1e-12
        )


def test_ctf_even_and_units():
    q = np.array([[0, 0], [0.01, 0.02], [0.05, 0.07]])
    p = np.array([[128, 2, 15000, 16000, 35, 300, 2.7, 0.1, 0]])
    np.testing.assert_allclose(ctf(q, p), ctf(-q, p), atol=1e-6)
    np.testing.assert_allclose(ctf(q, p)[0, 0], -0.1)


def test_fsc_identity_noise_and_censoring():
    rng = np.random.default_rng(8)
    v = rng.normal(size=(32, 32, 32))
    f = np.fft.fftshift(np.fft.fftn(v))
    c = fsc(f, f, 2)
    np.testing.assert_allclose(c[:, 2], 1, atol=1e-6)
    assert resolution(c)["status"] == "censored_at_sampled_limit"
    noise = np.fft.fftshift(np.fft.fftn(rng.normal(size=v.shape)))
    assert np.nanmean(np.abs(fsc(f, noise, 2)[3:, 2])) < 0.08


def test_volume_fft_transforms_axial_dimension_and_roundtrips():
    from fourier_splats.physics import fft_volume_center, volume_from_fourier
    n = 16
    z = np.arange(-n // 2, n // 2)
    v = np.broadcast_to(np.cos(2 * np.pi * 3 * z[:, None, None] / n), (n, n, n))
    f = fft_volume_center(v)
    nonzero = np.argwhere(np.abs(f) > 1e-8)
    np.testing.assert_array_equal(nonzero, [[n // 2 - 3, n // 2, n // 2], [n // 2 + 3, n // 2, n // 2]])
    np.testing.assert_allclose(volume_from_fourier(f), v, atol=1e-7)


def test_grid_inverse_real():
    rng = np.random.default_rng(8)
    a, b = design(np.zeros((1, 3)), 16)
    r = rng.normal(size=a.shape[1])
    im = rng.normal(size=a.shape[1])
    f = evaluate_grid(r, im, 16)
    v = np.fft.ifftn(np.fft.ifftshift(f))
    assert np.max(np.abs(v.imag)) < 1e-7


def test_adaptive_reference_and_geometry_gradients():
    import torch
    from fourier_splats.adaptive import FourierGaussianPairs

    rng = np.random.default_rng(77)
    centers = rng.normal(size=(4, 3))
    coeff = rng.normal(size=4) + 1j * rng.normal(size=4)
    L = np.tile(np.eye(3), (4, 1, 1))
    L[:, 1, 0] = 0.2
    model = FourierGaussianPairs(centers, coeff, L)
    k = torch.tensor(rng.normal(size=(10, 3)), dtype=torch.float64)
    got = model(k)
    want = hermitian_pairs(k.numpy(), centers, L @ L.transpose(0, 2, 1), coeff)
    np.testing.assert_allclose(
        got.detach().numpy()[:, 0], want.real, rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        got.detach().numpy()[:, 1], want.imag, rtol=1e-12, atol=1e-12
    )
    loss = (got**2).sum()
    loss.backward()
    analytic = model.centers.grad[0, 0].item()
    with torch.no_grad():
        original = model.centers[0, 0].item()
        eps = 1e-5
        model.centers[0, 0] = original + eps
        p = (model(k) ** 2).sum().item()
        model.centers[0, 0] = original - eps
        m = (model(k) ** 2).sum().item()
    np.testing.assert_allclose(analytic, (p - m) / (2 * eps), rtol=1e-6)
