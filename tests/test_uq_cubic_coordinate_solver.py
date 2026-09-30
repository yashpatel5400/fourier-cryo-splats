"""Check coordinate optimization against the independent conic formulation."""
import json
import numpy as np
import pytest
from test_uq_cubic_design_review import two_particle_fixture, conic_with_explicit_particle_packing
from fourier_splats.uq_block_preconditioner import BlockLowRankCoordinates
from fourier_splats.uq_cubic_coordinate_solver import optimize_with_coordinates


def nonsymmetric_coordinates(obj):
    rng = np.random.default_rng(650181)
    a = rng.normal(size=(obj.n, 2, obj.nq, obj.nq))
    blocks = a@a.swapaxes(-1, -2)+.2*np.eye(obj.nq)
    factors = [rng.normal(size=(obj.n*obj.nq, 3)) for _ in range(2)]
    return BlockLowRankCoordinates(blocks, factors)


def test_coordinate_solver_matches_two_particle_conic_optimum():
    obj, initial = two_particle_fixture()
    optimum, _, status = conic_with_explicit_particle_packing(obj, solver='SCS')
    assert status == 'optimal'
    transform = nonsymmetric_coordinates(obj)
    T = np.column_stack([transform.to_weights(e) for e in np.eye(initial.size)])
    assert np.linalg.norm(T-T.T) > .01
    result = optimize_with_coordinates(obj, initial, coordinates=transform,
        optimization_seed=650182, certificate_seed=650183, modes_count=1,
        ritz_subspace=7, ritz_tolerance=1e-9, maxiter=250, max_evaluations=300,
        ftol=1e-13, gtol=1e-9, power_iterations=80)
    assert result['optimizer_success']
    assert abs(result['selected_approximate_objective']-optimum) < 2e-5
    assert result['dual_lower_bound'] <= optimum+3e-7
    assert result['sum_objective_upper'] >= optimum-3e-7
    assert result['spectral_upper_certificate']['seed'] == 650183
    print('COORDINATE_CHECK', json.dumps({k: result[k] for k in [
        'selected_approximate_objective', 'dual_lower_bound', 'sum_objective_upper',
        'relative_sum_gap', 'optimizer_success', 'optimizer_iterations']}), 'CONIC', optimum)


@pytest.mark.parametrize('design, certificate', [(1, 1), (True, 2), (2, -1), (2, 1.5)])
def test_separate_explicit_streams_required_before_optimization(design, certificate):
    obj, initial = two_particle_fixture()
    with pytest.raises(ValueError):
        optimize_with_coordinates(obj, initial, coordinates=nonsymmetric_coordinates(obj),
            optimization_seed=design, certificate_seed=certificate)
