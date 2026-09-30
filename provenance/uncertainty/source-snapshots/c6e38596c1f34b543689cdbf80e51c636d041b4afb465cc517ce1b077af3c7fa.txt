import numpy as np
from scipy.linalg import block_diag
from scipy.spatial.transform import Rotation
from fourier_splats.uq_fourier_variational import HermitianTrilinearOperator
from fourier_splats.uq_fourier_pose_baseline import trilinear_pose_jacobian, GaussianNuisanceWhitening
from fourier_splats.uq_baselines import gaussian_reference_operator


def test_trilinear_rotation_translation_jacobian_against_nonlinear_forward():
    rng = np.random.default_rng(670101); n, nq = 3, 5
    q = rng.uniform(-.9, .9, (n, nq, 2))
    k = np.pad(q, ((0, 0), (0, 0), (0, 1)))@Rotation.random(n, random_state=rng).as_matrix()
    ctf = rng.normal(size=(n, nq)); noise = .3
    op = HermitianTrilinearOperator(k, ctf, noise, box=5)
    x = rng.normal(size=op.shape[1]); J, _ = trilinear_pose_jacobian(op, k, q, ctf, noise, x)
    directions = rng.normal(size=(n, 5)); step = 1e-6
    def forward(t):
        rotated = k@Rotation.from_rotvec(t*directions[:, :3]).as_matrix()
        value = (HermitianTrilinearOperator(rotated, ctf, noise, box=5).matrix@x).reshape(n, 2*nq)
        z = (value[:, :nq]+1j*value[:, nq:])*np.exp(-2j*np.pi*t*np.einsum('nqi,ni->nq', q, directions[:, 3:]))
        return np.concatenate([z.real, z.imag], axis=1)
    np.testing.assert_allclose((forward(step)-forward(-step))/(2*step),
        np.einsum('nmp,np->nm', J, directions), rtol=3e-8, atol=2e-8)


def test_marginal_whitening_matches_dense_joint_gaussian_posterior():
    rng = np.random.default_rng(670102); n, m, p, d = 3, 7, 5, 2
    j = rng.normal(size=(n, m, d)); a = rng.normal(size=(n*m, p))
    fullj = block_diag(*j); tau = .8; ell = rng.normal(size=p)
    transform = GaussianNuisanceWhitening(j)
    covariance = np.eye(n*m)+fullj@fullj.T
    W = np.column_stack([transform.apply(e) for e in np.eye(n*m)])
    np.testing.assert_allclose(W@covariance@W, np.eye(n*m), atol=3e-14)
    np.testing.assert_allclose(W, W.T, atol=1e-15)
    fit = gaussian_reference_operator(transform.whiten_operator(a), ell, tau,
        gram_diagonal=np.sum(a*a, axis=0), rtol=1e-12)
    weights = transform.apply(fit['weights'])
    joint = np.column_stack([a, fullj])
    precision = joint.T@joint+np.diag([1/tau**2]*p+[1.]*(n*d))
    target = np.r_[ell, np.zeros(n*d)]
    independent = np.linalg.solve(precision, target)
    np.testing.assert_allclose(weights, joint@independent, atol=2e-13)
    error = a.T@weights-ell
    variance = weights@covariance@weights+tau**2*(error@error)
    np.testing.assert_allclose(variance, target@independent, atol=2e-13)
    fixed_variance = ell@np.linalg.solve(a.T@a+np.eye(p)/tau**2, ell)
    assert variance >= fixed_variance
    zero = GaussianNuisanceWhitening(np.zeros_like(j))
    np.testing.assert_array_equal(zero.apply(np.arange(n*m)), np.arange(n*m))
