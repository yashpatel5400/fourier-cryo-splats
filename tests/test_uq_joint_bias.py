import numpy as np
from numpy.polynomial.legendre import leggauss
from fourier_splats.uq_joint_bias import gaussian_fourier_moments,target_pose_cross,joint_density_pose_bias,MONOMIALS,sharp_cube_cubic_coefficients
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator


def test_gaussian_moments_match_independent_quadrature_with_boundary_terms():
    rng=np.random.default_rng(610121);k=rng.normal(size=(7,3));centers=[[.32,-.12,.2],[-.1,.27,0]];signs=[1.,-.7];width=.19
    nodes,weights=leggauss(55);nodes/=2;weights/=2
    expected=np.zeros((len(k),10),complex)
    for mu,sign in zip(centers,signs):
        for j,powers in enumerate(MONOMIALS):
            product=np.ones(len(k),complex)
            for axis,power in enumerate(powers):
                density=np.exp(-(nodes-mu[axis])**2/(2*width**2))/(np.sqrt(2*np.pi)*width)
                product*=np.exp(2j*np.pi*k[:,axis,None]*nodes)@(weights*density*nodes**power)
            expected[:,j]+=sign*product
    np.testing.assert_allclose(gaussian_fourier_moments(k,centers,signs,width),expected,atol=2e-13,rtol=2e-12)


def test_analytic_target_pose_cross_matches_direct_field_quadrature():
    rng=np.random.default_rng(610122);k=.7*rng.normal(size=(2,3,3));q=rng.normal(size=(2,3,2));ctf=rng.normal(size=(2,3));w=rng.normal(size=12)
    op=PolynomialPoseFieldOperator(k,q,ctf,w,.8,.03,.002,order=28,backend='direct')
    centers=[[.1,-.1,0]];width=.16
    ell=np.exp(-np.sum((op.xyz-centers[0])**2,axis=1)/(2*width**2))/(2*np.pi*width**2)**1.5
    np.testing.assert_allclose(target_pose_cross(op,centers,[1],width),op.rmatvec(op.sqrt_quad*ell),atol=2e-12,rtol=2e-10)


def test_joint_bias_bounds_actual_hilbert_pose_error_and_improves_orthogonal_case():
    rng=np.random.default_rng(610123);F=rng.normal(size=(8,4));h=rng.normal(size=8);pilot=rng.normal(size=8)
    L=1.7;B=2.;P=np.linalg.norm(pilot);upper=L*np.linalg.norm(F,2)
    bound=joint_density_pose_bias(np.linalg.norm(h),upper,np.linalg.norm(F.T@h),L,B,P,0)
    for _ in range(200):
        v=rng.normal(size=4);v*=L/np.linalg.norm(v)
        assert abs(pilot@(F@v))+B*np.linalg.norm(h-F@v)<=bound['bias_upper']+1e-10
    orthogonal=joint_density_pose_bias(1.,1.,0.,1.,2.,1.,0.)
    np.testing.assert_allclose(orthogonal['bias_upper'],2*np.sqrt(2)+1)
    assert orthogonal['bias_upper']<orthogonal['triangle_bias_upper']


def test_sharp_cube_moments_by_independent_polynomial_integration():
    rng=np.random.default_rng(610141);nodes,weights=leggauss(5);nodes/=2;weights/=2
    xyz=np.stack(np.meshgrid(nodes,nodes,nodes,indexing='ij'),axis=-1).reshape(-1,3)
    qw=np.prod(np.stack(np.meshgrid(weights,weights,weights,indexing='ij'),axis=-1),axis=-1).ravel()
    u=rng.normal(size=(200,3));u/=np.linalg.norm(u,axis=1)[:,None];u=np.vstack([u,np.ones(3)/np.sqrt(3)])
    fourth=(xyz@u.T)**4; sixth=(xyz@u.T)**6
    assert np.max(qw@fourth)<=13/720+1e-14
    assert np.max(qw@sixth)<=205/36288+1e-14
    np.testing.assert_allclose(qw@fourth[:,-1],13/720,rtol=1e-12)
    np.testing.assert_allclose(qw@sixth[:,-1],205/36288,rtol=1e-12)


def test_sharp_cubic_integrated_bound_covers_exact_single_frequency_remainder():
    from scipy.spatial.transform import Rotation
    from fourier_splats.uq_continuous_pose import pose_derivative_fields,cube_quadrature
    from fourier_splats.uq_pose_optimization import cubic_penalty_coefficients
    rng=np.random.default_rng(610142);xyz,qw=cube_quadrature(24)
    for _ in range(15):
        k=rng.normal(size=(1,3));q=rng.normal(size=(1,2));angle=.08;shift=.01
        u=rng.normal(size=5);u/=np.linalg.norm(u);coefficient=np.exp(1j*rng.uniform(-np.pi,np.pi))
        first,second=pose_derivative_fields(k,q,np.array([coefficient]),xyz,angle,shift,backend='direct')
        from fourier_splats.uq_continuous_pose import PAIRS,PAIR_SCALE
        polynomial=first@u+.5*second@np.array([u[a]*u[b]*s for (a,b),s in zip(PAIRS,PAIR_SCALE)])
        rotated=k@Rotation.from_rotvec(angle*u[:3]).as_matrix()
        exact=(coefficient*(np.exp(2j*np.pi*(xyz@rotated.T).ravel()+2j*np.pi*shift*(q@u[3:]).item())-np.exp(2j*np.pi*(xyz@k.T).ravel()))).real
        remainder=np.sqrt(qw@(exact-polynomial)**2)
        sharp=sharp_cube_cubic_coefficients(k[None],q[None],np.ones((1,1)),angle,shift,1.)[0,0]
        old=cubic_penalty_coefficients(k[None],q[None],np.ones((1,1)),angle,shift,1.)[0,0]
        assert remainder<=sharp+1e-11
        assert sharp<old
