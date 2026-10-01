import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp
from fourier_splats.uq_pose_importance import PoseProposal,importance_summary


def test_acg_full_covariance_and_antipodes():
    rng=np.random.default_rng(624);centers=Rotation.random(2,random_state=rng).as_quat()
    cov=np.array([np.diag([.1,.3,.7]),np.diag([.2,.8,.4])]);weights=np.array([.3,.7])
    model=PoseProposal(centers,cov,weights);points=Rotation.random(13,random_state=rng).as_quat()
    density=np.full(13,.05)
    basis=np.eye(4)
    for i in range(2):
        # Independent full four-dimensional covariance in the global quaternion frame.
        transform=(Rotation.from_quat(centers[i])*Rotation.from_quat(basis)).as_quat().T
        local=np.eye(4);local[:3,:3]=cov[i]/4
        full=transform@local@transform.T
        quadratic=np.einsum('ni,ij,nj->n',points,np.linalg.inv(full),points)
        density+=.95*weights[i]*np.linalg.det(full)**(-.5)*quadratic**(-2)
    np.testing.assert_allclose(model.log_density(points),np.log(density),atol=2e-14)
    np.testing.assert_allclose(model.log_density(-points),model.log_density(points),atol=2e-14)


def test_uniform_and_inverse_density_sampling():
    center=Rotation.from_rotvec([.4,-.5,.7]).as_quat()[None]
    isotropic=PoseProposal(center,np.eye(3)[None]*4,np.ones(1))
    rng=np.random.default_rng(382);points,_=isotropic.sample(1024,rng)
    np.testing.assert_allclose(isotropic.log_density(points),0,atol=2e-15)
    model=PoseProposal(center,np.diag([.4,.7,.9])[None],np.ones(1),haar_weight=.2)
    points,_=model.sample(100000,rng);z=np.exp(-model.log_density(points))
    assert abs(z.mean()-1)<6*z.std(ddof=1)/np.sqrt(len(z))
    f=points[:,3]**2*z
    assert abs(f.mean()-.25)<6*f.std(ddof=1)/np.sqrt(len(f))


def test_importance_log_integrals_against_direct_weights():
    k=np.array([[-10.,-11.],[-12.,-10.],[-8.,-9.]]);q=np.log([.5,1.,2.])
    r=importance_summary(k,q);w=np.exp(k)/np.exp(q[:,None]);p=w/w.sum(0)
    np.testing.assert_allclose(r['log_integrals'],np.log(w.mean(0)),atol=2e-15)
    np.testing.assert_allclose(r['ess'],1/np.sum(p*p,axis=0),atol=2e-15)
    np.testing.assert_allclose(r['log_ratio'],np.log(w[:,1].sum()/w[:,0].sum()),atol=2e-15)
