import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_cell_pose_pair import CellPairDistance, refine_cell_density
from fourier_splats.uq_continuous import cell_forward


def test_pose_pair_forward_and_gradient_against_independent_direct_cell_integrals():
    rng = np.random.default_rng(610421); box = 4; n, nq = 2, 3
    k = rng.normal(size=(n, nq, 3)); q = rng.normal(size=(n, nq, 2)); ctf = rng.normal(size=(n, nq))
    densities = rng.normal(size=(2, box**3)); angle, shift, noise = .2, .08, .7
    problem = CellPairDistance(k, q, ctf, densities, box, noise, angle, shift)
    xyz = (np.indices((box,)*3).reshape(3, -1).T[:, ::-1]+.5)/box-.5
    def direct(coordinates):
        u = problem.project(coordinates); predictions = []
        for j in range(2):
            rotated = k@Rotation.from_rotvec(angle*u[j, :, :3]).as_matrix()
            phase = shift*np.einsum('nqa,na->nq', q, u[j, :, 3:])
            f = np.exp(-2j*np.pi*(rotated@xyz.T))@densities[j]/box**1.5
            predictions.append(ctf/noise*f*np.prod(np.sinc(rotated/box), axis=-1)*np.exp(-2j*np.pi*phase))
        return float(np.sum(abs(predictions[1]-predictions[0])**2))
    for scale in [.15, 1.2]:
        coordinates = rng.normal(size=(2, n, 5))*scale
        value, gradient = problem.value_gradient(coordinates)
        np.testing.assert_allclose(value, direct(coordinates), rtol=3e-12, atol=1e-12)
        flat = coordinates.ravel(); numerical = []
        for a in range(len(flat)):
            direction = np.eye(len(flat))[a]*1e-5
            numerical.append((direct(flat+direction)-direct(flat-direction))/(2e-5))
        np.testing.assert_allclose(gradient, numerical, rtol=3e-6, atol=2e-8)
        distance, check = problem.independent_distance(coordinates)
        np.testing.assert_allclose(check['unpadded_distance']**2, value, rtol=3e-12)
        assert distance >= np.sqrt(direct(coordinates))
        assert np.linalg.norm(problem.project(coordinates), axis=-1).max() <= 1+1e-14


def test_nested_pilot_refinement_preserves_norm_and_exact_forward_field():
    rng = np.random.default_rng(610423); coarse = rng.normal(size=4**3)
    fine = refine_cell_density(coarse, 4, 12)
    np.testing.assert_allclose(np.linalg.norm(fine), np.linalg.norm(coarse), rtol=1e-15)
    k = rng.normal(size=(2, 5, 3)); ctf = rng.normal(size=(2, 5))
    np.testing.assert_allclose(cell_forward(k, ctf, fine, 12, .4), cell_forward(k, ctf, coarse, 4, .4), rtol=2e-9, atol=1e-10)
