import numpy as np
from scipy.fft import fftn, ifftn, fftshift, ifftshift
from fourier_splats.uq_data import VoxelReference,VoxelObservationOperator,discrete_density_functionals,local_weights
from fourier_splats.uq_physics import gaussian_pair_features


def test_nufft_generator_matches_direct_sum_and_cartesian_fft():
    rng=np.random.default_rng(204);box=8
    v=rng.normal(size=(box,)*3);ref=VoxelReference(v,field_A=40)
    q=np.array([[0.,0.,0.],[1.,2.,-1.],[-2.,1.,0.],[.21,-1.43,.78]])
    actual=ref.fourier(q)
    g=np.arange(-box//2,box//2);z,y,x=np.meshgrid(g,g,g,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    expected=np.exp(-2j*np.pi*q@xyz.T/box)@v.ravel()
    assert np.allclose(actual,expected,rtol=1e-8,atol=1e-8)
    f=fftshift(fftn(ifftshift(v)))
    assert np.allclose(actual[:3],f[q[:3,2].astype(int)+4,q[:3,1].astype(int)+4,q[:3,0].astype(int)+4])


def test_functionals_match_rendered_gaussian_map():
    rng=np.random.default_rng(205);box=12
    centers=rng.normal(size=(6,3));c=rng.normal(size=12)
    weights=np.stack([local_weights(box,[1,2,-1],1.4),local_weights(box,[-2,-1,2],2)])
    ell=discrete_density_functionals(weights,centers,.8)
    g=np.arange(-box//2,box//2);z,y,x=np.meshgrid(g,g,g,indexing='ij')
    k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    f=(gaussian_pair_features(k,centers,.8)@c).reshape((box,)*3)
    rendered=fftshift(ifftn(ifftshift(f))).real
    expected=np.sum(weights*rendered,axis=(1,2,3))
    assert np.allclose(ell@c,expected,atol=1e-12)


def test_full_grid_adjoint_and_outside_dictionary_bias():
    rng=np.random.default_rng(901);n,q,box=3,7,6
    k=rng.normal(size=(n,q,3));transfer=rng.normal(size=(n,q))
    op=VoxelObservationOperator(k,transfer,box,.4)
    v=rng.normal(size=box**3);w=rng.normal(size=n*2*q)
    assert np.allclose(w@op.forward(v),v@op.adjoint(w),rtol=1e-10,atol=1e-9)
    grid=np.arange(-box//2,box//2);z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
    xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    complex_design=np.exp(-2j*np.pi*k@xyz.T/box)*transfer[...,None]/.4
    design=np.concatenate([complex_design.real,complex_design.imag],axis=1).reshape(n*2*q,-1)
    assert np.allclose(op.forward(v),design@v,rtol=1e-9,atol=1e-8)
    assert np.allclose(op.adjoint(w),design.T@w,rtol=1e-9,atol=1e-8)
    # Exact ambient support function remains valid for a perturbation excluded
    # from a low-rank density model; its projection need not bound this bias.
    ell=rng.normal(size=box**3);basis=np.linalg.qr(rng.normal(size=(box**3,8)))[0]
    residual=op.bias_residual(w,ell);outside=residual-basis@(basis.T@residual)
    delta=.7*outside/np.linalg.norm(outside)
    bias=w@op.forward(delta)-ell@delta
    assert abs(bias)<=.7*np.linalg.norm(residual)+1e-8
    assert abs(bias)>.7*np.linalg.norm(basis.T@residual)


def test_ambient_pose_terms_against_nonlinear_voxel_projections():
    from fourier_splats.uq_data import SupportedVoxelOperator,voxel_pose_terms
    from scipy.spatial.transform import Rotation
    rng=np.random.default_rng(926);n,q,box=3,7,6
    k=rng.normal(size=(n,q,3));detector=rng.normal(size=(n,q,2));ctf=rng.normal(size=(n,q))
    mask=rng.uniform(size=box**3)>.4
    base=VoxelObservationOperator(k,ctf,box,.7);op=SupportedVoxelOperator(base,mask)
    pilot=rng.normal(size=op.shape[1]);pilot/=np.linalg.norm(pilot)
    B=.4;angle=.08;shift=.15
    terms=voxel_pose_terms(op,pilot,detector,angle,shift,B)
    def forward(coef,u):
        rotations=Rotation.from_rotvec(angle*u[:,:3]).as_matrix()
        kp=np.einsum('nqi,nij->nqj',k,rotations)
        shifted=VoxelObservationOperator(kp,ctf,box,.7).forward(op.expand(coef)).reshape(n,2*q)
        z=shifted[:,:q]+1j*shifted[:,q:]
        z*=np.exp(-2j*np.pi*shift*np.einsum('nqi,ni->nq',detector,u[:,3:])/box)
        return np.concatenate([z.real,z.imag],axis=1)
    for axis in range(5):
        u=np.zeros((n,5));u[:,axis]=1e-5
        numeric=(forward(pilot,u)-forward(pilot,-u))/(2e-5)
        assert np.allclose(numeric,terms['jacobian'][:,:,axis],rtol=1e-6,atol=1e-8)
    for _ in range(20):
        delta=rng.normal(size=len(pilot));delta*=B/np.linalg.norm(delta)
        u=rng.normal(size=(n,5));u/=np.linalg.norm(u,axis=1,keepdims=True)
        delta_terms=voxel_pose_terms(op,delta,detector,angle,shift,B)
        linear=np.einsum('nmq,nq->nm',terms['jacobian']+delta_terms['jacobian'],u)
        remainder=forward(pilot+delta,u)-op.forward(pilot+delta).reshape(n,2*q)-linear
        assert np.all(np.linalg.norm(remainder,axis=1)<=terms['remainder']+1e-8)
        w=rng.normal(size=(n,2*q))
        interaction=np.einsum('nmq,nq,nm->n',delta_terms['jacobian'],u,w)
        bound=B*np.linalg.norm(np.einsum('nmq,nm->nq',terms['interaction_factor'],w),axis=1)
        assert np.all(np.abs(interaction)<=bound+1e-8)
    design=op.forward_columns(np.eye(op.shape[1]))
    assert np.allclose(op.gram_diagonal,np.sum(design*design,axis=0),rtol=1e-9)


def test_ambient_quadratic_pose_derivatives_and_uniform_cubic_bound():
    from fourier_splats.uq_data import SupportedVoxelOperator,voxel_pose_terms
    from scipy.spatial.transform import Rotation
    rng=np.random.default_rng(930);n,q,box=2,6,6
    k=rng.normal(size=(n,q,3));detector=rng.normal(size=(n,q,2));ctf=rng.normal(size=(n,q))
    op=SupportedVoxelOperator(VoxelObservationOperator(k,ctf,box,.7),rng.uniform(size=box**3)>.5)
    pilot=rng.normal(size=op.shape[1]);pilot/=np.linalg.norm(pilot)
    B=.4;angle=.16;shift=.2
    terms=voxel_pose_terms(op,pilot,detector,angle,shift,B,order=2)
    def forward(coef,u):
        kp=np.einsum('nqi,nij->nqj',k,Rotation.from_rotvec(angle*u[:,:3]).as_matrix())
        data=VoxelObservationOperator(kp,ctf,box,.7).forward(op.expand(coef)).reshape(n,2*q)
        z=data[:,:q]+1j*data[:,q:];z*=np.exp(-2j*np.pi*shift*np.einsum('nqi,ni->nq',detector,u[:,3:])/box)
        return np.concatenate([z.real,z.imag],axis=1)
    eps=.002
    for a,b in [(0,0),(0,1),(1,2),(0,3),(3,4),(4,4)]:
        ua=np.zeros((n,5));ub=np.zeros((n,5));ua[:,a]=eps;ub[:,b]=eps
        numeric=(forward(pilot,ua+ub)-forward(pilot,ua-ub)-forward(pilot,-ua+ub)+forward(pilot,-ua-ub))/(4*eps**2)
        np.testing.assert_allclose(numeric,terms['quadratic_jacobian'][:,:,a,b],rtol=2e-4,atol=2e-7)
    for _ in range(15):
        delta=rng.normal(size=len(pilot));delta*=B/np.linalg.norm(delta)
        u=rng.normal(size=(n,5));u/=np.linalg.norm(u,axis=1,keepdims=True)
        pert=voxel_pose_terms(op,delta,detector,angle,shift,B,order=2)
        linear=np.einsum('nma,na->nm',terms['jacobian']+pert['jacobian'],u)
        quadratic=.5*np.einsum('nmab,na,nb->nm',terms['quadratic_jacobian']+pert['quadratic_jacobian'],u,u)
        remainder=forward(pilot+delta,u)-op.forward(pilot+delta).reshape(n,2*q)-linear-quadratic
        assert np.all(np.linalg.norm(remainder,axis=1)<=terms['remainder']+1e-8)
        w=rng.normal(size=(n,2*q));interaction=.5*np.einsum('nm,nmab,na,nb->n',w,pert['quadratic_jacobian'],u,u)
        bound=.5*B*np.linalg.norm(np.einsum('nmk,nm->nk',terms['quadratic_interaction_factor'],w),axis=1)
        assert np.all(np.abs(interaction)<=bound+1e-8)
