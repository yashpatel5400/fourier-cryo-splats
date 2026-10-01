import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_pose_importance import PoseProposal
from fourier_splats.uq_pose_catalog import CatalogPoseProposal


def model():
    rng=np.random.default_rng(328)
    local=PoseProposal(Rotation.random(2,random_state=rng).as_quat(),
                       np.array([np.diag([.2,.4,.6]),np.diag([.3,.8,.5])]),[.4,.6])
    centers=Rotation.random(3,random_state=rng).as_quat()
    return CatalogPoseProposal(local,centers,[.2,.3,.5],.7,masses=(.15,.35,.5))


def test_catalog_density_independent_global_covariance():
    m=model();rng=np.random.default_rng(625)
    points=Rotation.random(17,random_state=rng).as_quat()
    density=np.full(len(points),m.masses[0])
    covs=list(m.local.covariances/4)+[np.eye(3)*m.s**2]*len(m.centers)
    centers=list(m.local.centers.as_quat())+list(m.centers)
    weights=np.r_[m.masses[1]*m.local.weights,m.masses[2]*m.weights]
    for center,cov,weight in zip(centers,covs,weights):
        transform=(Rotation.from_quat(center)*Rotation.from_quat(np.eye(4))).as_quat().T
        local=np.eye(4);local[:3,:3]=cov
        global_cov=transform@local@transform.T
        z=np.einsum('ni,ij,nj->n',points,np.linalg.inv(global_cov),points)
        density+=weight*np.linalg.det(global_cov)**(-.5)*z**(-2)
    np.testing.assert_allclose(m.log_density(points,chunk=3),np.log(density),atol=3e-14)
    np.testing.assert_allclose(m.log_density(-points),m.log_density(points),atol=3e-14)


def test_catalog_sampling_against_haar_moments():
    m=model();points,family,component=m.sample(100000,np.random.default_rng(729))
    w=np.exp(-m.log_density(points))
    for f,target in [(np.ones(len(w)),1),(points[:,0]**2,.25),(points[:,2]**4,.125)]:
        z=f*w
        assert abs(z.mean()-target)<6*z.std(ddof=1)/np.sqrt(len(z))
    np.testing.assert_allclose(np.bincount(family,minlength=3)/len(family),m.masses,atol=.005)
    assert np.all(component[family==0]==-1)
    one,_,_=m.sample(1,np.random.default_rng(6))
    assert one.shape==(1,4)


def test_uniform_catalog_limit_and_small_kernel_center():
    center=np.array([[0.,0.,0.,1.]])
    local=PoseProposal(center,np.eye(3)[None]*4,[1.])
    points=Rotation.random(100,random_state=654).as_quat()
    m=CatalogPoseProposal(local,center,[1.],2.)
    np.testing.assert_allclose(m.log_density(points),0,atol=3e-15)
    m=CatalogPoseProposal(local,center,[1.],.01)
    # The quaternion kernel at its centre has exactly the determinant factor.
    expected=m.masses[0]+m.masses[1]+m.masses[2]*m.s**-3
    np.testing.assert_allclose(m.log_density(center),np.log(expected),atol=3e-14)
