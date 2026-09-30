import numpy as np
from numpy.polynomial.legendre import leggauss
from fourier_splats.uq_joint_bias import gaussian_fourier_moments,target_pose_cross,joint_density_pose_bias,MONOMIALS,sharp_cube_cubic_coefficients,residual_pose_cross_bound,scaled_pose_radius,pose_scale_record
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


def test_archive_gaussian_moments_against_200_node_separable_rule():
    k=np.array([[5,0,0],[-5,0,0],[0,3,4],[5/np.sqrt(3)]*3,[2,-3,np.sqrt(12)]])
    width=.07;centers=[[0,0,.08],[0,0,-.08]];signs=[1.,-1.]
    x,weight=leggauss(200);x/=2;weight/=2;expected=np.zeros((len(k),10),complex)
    for mu,sign in zip(centers,signs):
        for j,powers in enumerate(MONOMIALS):
            integral=np.ones(len(k),complex)
            for axis,power in enumerate(powers):
                density=np.exp(-(x-mu[axis])**2/(2*width**2))/(np.sqrt(2*np.pi)*width)
                integral*=np.exp(2j*np.pi*k[:,axis,None]*x)@(weight*density*x**power)
            expected[:,j]+=sign*integral
    np.testing.assert_allclose(gaussian_fourier_moments(k,centers,signs,width),expected,rtol=2e-11,atol=1e-14)


def test_archive_cross_quadrature_pad_against_order96_direct_evaluation():
    from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
    k=np.array([[[5.,0,0],[3.,4,0]]]);q=k[:,:,:2].copy();ctf=np.array([[.8,-.6]])
    w=np.array([.7,-.9,.3,.4]);noise=.8;angle=np.deg2rad(2);shift=.5/236.16
    low=PolynomialPoseFieldOperator(k,q,ctf,w,noise,angle,shift,order=32,backend='direct')
    bound=residual_pose_cross_bound(low,w,noise,[[0,0,0]],[1],.07)
    high=PolynomialPoseFieldOperator(k,q,ctf,w,noise,angle,shift,order=96,backend='direct')
    # Assemble the field by explicit exponentials, independently of NUFFT.
    c=(w[:2]+1j*w[2:])*ctf[0]/noise
    field=(np.exp(2j*np.pi*(high.xyz@k[0].T))@c).real
    exact=target_pose_cross(high,[[0,0,0]],[1],.07)-high.rmatvec(high.sqrt_quad*field)
    error=np.linalg.norm(bound['cross_vector']-exact)
    assert error<=bound['quadrature_norm_pad']
    assert np.linalg.norm(exact)<=bound['norm_upper']


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


def test_archive_cubic_remainder_at_maximizing_diagonal_direction():
    from scipy.spatial.transform import Rotation
    from fourier_splats.uq_continuous_pose import pose_derivative_fields,cube_quadrature,PAIRS,PAIR_SCALE
    direction=np.ones(3)/np.sqrt(3);kh=np.array([1.,-1.,0.])/np.sqrt(2)
    omega=np.cross(direction,kh);np.testing.assert_allclose(np.cross(kh,omega),direction,atol=1e-15)
    k=5*kh[None];q=np.array([[5.,0]]);angle=np.deg2rad(2);shift=.5/236.16
    u=np.r_[omega,0.,0.];xyz,qw=cube_quadrature(48);c=np.array([np.exp(.7j)])
    first,second=pose_derivative_fields(k,q,c,xyz,angle,shift,backend='direct')
    polynomial=first@u+.5*second@np.array([u[a]*u[b]*s for (a,b),s in zip(PAIRS,PAIR_SCALE)])
    rotated=k@Rotation.from_rotvec(angle*omega).as_matrix()
    exact=(c*(np.exp(2j*np.pi*(xyz@rotated.T))-np.exp(2j*np.pi*(xyz@k.T)))).real.ravel()
    remainder=np.sqrt(qw@(exact-polynomial)**2)
    bound=sharp_cube_cubic_coefficients(k[None],q[None],np.ones((1,1)),angle,shift,1.)[0,0]
    assert remainder<=bound


def test_joint_and_product_pose_lifts_require_different_radii():
    from fourier_splats.uq_continuous_pose import PAIRS,PAIR_SCALE
    d=np.array([.8,.3,1.1,.6]);xi=np.array([[1.,0,0,1.,0],[0,1.,0,0,1.]])
    lift=np.concatenate([(np.sqrt(d[:2,None])*xi).ravel(),
        (np.sqrt(d[2:,None])*np.array([[row[a]*row[b]*s for (a,b),s in zip(PAIRS,PAIR_SCALE)] for row in xi])).ravel()])
    assert np.linalg.norm(lift)>scaled_pose_radius(d,'joint')
    np.testing.assert_allclose(np.linalg.norm(lift),scaled_pose_radius(d,'product'),rtol=1e-14)
    xi/=np.sqrt(2)
    joint_norm2=np.sum(d[:2]*np.linalg.norm(xi,axis=1)**2)+np.sum(d[2:]*np.linalg.norm(xi,axis=1)**4)
    np.testing.assert_allclose(joint_norm2,scaled_pose_radius(d,'joint')**2,rtol=1e-14)
    denom=np.sqrt(np.r_[np.repeat(d[:2],5),np.repeat(d[2:],15)])
    saved=pose_scale_record(d,denom)
    # Same sum is insufficient: permuting blocks changes the operator metadata.
    changed=pose_scale_record(d[::-1],denom[::-1])
    assert saved['group_scales_sha256']!=changed['group_scales_sha256']
