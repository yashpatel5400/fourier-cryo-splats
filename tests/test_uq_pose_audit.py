import numpy as np
from fourier_splats.uq_data import VoxelObservationOperator,SupportedVoxelOperator,voxel_pose_terms
from fourier_splats.uq_pose_audit import pose_derivative_adjoints,audit_pose_weights


def test_moment_adjoints_match_explicit_column_derivatives_and_bounds():
    rng=np.random.default_rng(609340);n,q,box=3,9,6
    k=rng.normal(size=(n,q,3));detector=rng.normal(size=(n,q,2));ctf=rng.normal(size=(n,q))
    op=SupportedVoxelOperator(VoxelObservationOperator(k,ctf,box,.6),rng.uniform(size=box**3)>.4)
    pilot=rng.normal(size=op.shape[1]);pilot/=np.linalg.norm(pilot);w=rng.normal(size=(n,2*q))
    terms=voxel_pose_terms(op,pilot,detector,.13,.2,.5,order=2)
    for i in range(n):
        d1,d2=pose_derivative_adjoints(op,detector,w,.13,.2,i)
        fd1,fd2=pose_derivative_adjoints(op,detector,w,.13,.2,i,backend='nufft')
        np.testing.assert_allclose(fd1,d1,rtol=1e-8,atol=1e-9)
        np.testing.assert_allclose(fd2,d2,rtol=1e-8,atol=1e-9)
        np.testing.assert_allclose(pilot@d1,w[i]@terms['jacobian'][i],rtol=1e-9,atol=1e-9)
        np.testing.assert_allclose(np.einsum('p,pab->ab',pilot,d2),np.einsum('m,mab->ab',w[i],terms['quadratic_jacobian'][i]),rtol=1e-9,atol=1e-9)
        np.testing.assert_allclose(np.linalg.norm(d1),np.linalg.norm(w[i]@terms['interaction_factor'][i]),rtol=1e-9)
        np.testing.assert_allclose(np.linalg.norm(d2),np.linalg.norm(w[i]@terms['quadratic_interaction_factor'][i]),rtol=1e-9)
    audit=audit_pose_weights(op,pilot,detector,w,.13,.2,.5)
    np.testing.assert_allclose(audit['per_particle'][:,-1],terms['remainder']*np.linalg.norm(w,axis=1),rtol=1e-12)
