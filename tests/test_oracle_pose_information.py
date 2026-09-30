"""Verify physical units and derivatives against a separate nonlinear operator."""
import importlib.util
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_continuous_pose import pose_cell_forward

spec=importlib.util.spec_from_file_location('oracle_pose_information',Path(__file__).resolve().parents[1]/'scripts/audit_uq_oracle_pose_information.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_exact_cell_pose_derivatives_in_degrees_and_angstroms():
    rng=np.random.default_rng(711101);q=rng.normal(size=(2,7,2))
    k=np.pad(q,((0,0),(0,0),(0,1)))@Rotation.from_rotvec([.3,-.2,.4]).as_matrix()
    ctf=rng.normal(size=(2,7));rho=rng.normal(size=6**3);rho/=np.linalg.norm(rho)
    noise=.21;field=219.
    actual=module.cell_pose_jacobian(k,q,ctf,rho,6,noise,field)
    for a in range(5):
        delta=np.zeros((2,5));delta[:,a]=1e-4
        plus=pose_cell_forward(k,q,ctf,rho,6,noise,delta,np.pi/180,1/field).reshape(2,14)
        minus=pose_cell_forward(k,q,ctf,rho,6,noise,-delta,np.pi/180,1/field).reshape(2,14)
        np.testing.assert_allclose(actual[:,:,a],(plus-minus)/(2e-4),rtol=2e-6,atol=2e-9)


def test_information_scales_and_deficient_direction_is_not_hidden():
    rng=np.random.default_rng(711102);j=rng.normal(size=(1,18,5))
    first=module.local_information(j)[0];second=module.local_information(2*j)[0]
    np.testing.assert_allclose(second['covariance'],np.asarray(first['covariance'])/4,rtol=1e-12)
    assert first['known_shift_rotation_rms_degrees']<=first['rotation_rms_degrees']
    j[:,:,4]=j[:,:,3]
    bad=module.local_information(j)[0]
    assert bad['rank']==4 and bad['covariance'] is None
