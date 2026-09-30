import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from fourier_splats.uq_map_alignment import align_map_to_pilot, apply_map_alignment


@pytest.mark.parametrize('hand',[-1.,1.])
def test_alignment_on_independently_rendered_transformed_gaussians(hand):
    # Both volumes are analytic renders, so this does not validate interpolation
    # by generating the input with the same interpolation being tested.
    box = 24
    xyz = np.indices((box,)*3).reshape(3,-1).T-box/2
    centers = np.array([[-4.,-1.,0.],[3.,2.,-2.],[0.,-4.,4.],[3.,3.,4.]])
    amplitudes = [1.,.7,.4,.23]; widths = [1.6,1.9,1.3,1.1]
    rotation = Rotation.from_rotvec([.37,-.29,.44]).as_matrix()@np.diag([hand,1.,1.])
    shift = np.array([1.2,-.9,.7])
    def render(points):
        return sum(a*np.exp(-np.sum((points-c)**2,axis=1)/(2*w*w))
                   for a,c,w in zip(amplitudes,centers,widths)).reshape((box,)*3)
    pilot = render(xyz); source = render((xyz-shift)@rotation.T)
    result = align_map_to_pilot(source, pilot, 1., lowpass_A=3., box=16,max_evaluations=240)
    assert len(result['starts']) == 48
    assert sorted(x['hand_determinant'] for x in result['starts']) == [-1]*24+[1]*24
    transformed = apply_map_alignment(source,result['matrix'],result['offset'])
    assert np.corrcoef(pilot.ravel(),transformed.ravel())[0,1] > .985
    assert hand*np.linalg.det(result['matrix']) > 0
