"""Analytic fixed-pose audits in L2([-1/2,1/2]^3), beyond voxel subspaces.

The adjoint is a finite real Fourier sum. Its squared continuous norm follows
from the cube's sinc kernel, and Gaussian functional integrals use complex
normal CDFs. These formulas concern a fixed CTF, known pose and scalar whitening.
They do not extend the finite-voxel nonlinear pose guarantee by assertion.
"""
import numpy as np
from scipy.special import ndtr
from scipy.sparse.linalg import LinearOperator, cg
from scipy.stats import norm
from .uq_data import VoxelObservationOperator
from .uncertainty import bias_aware_half_width


def gaussian_target_integrals(k,centers,signs,width):
    """Fourier transform and squared L2 norm of a Gaussian contrast on the cube."""
    k=np.asarray(k);centers=np.asarray(centers);signs=np.asarray(signs);sigma=float(width)
    if sigma<=0:raise ValueError('Positive target width required')
    fourier=np.zeros(k.shape[:-1],complex)
    for mu,sign in zip(centers,signs):
        imaginary=2j*np.pi*sigma*k
        mass=ndtr((.5-mu)/sigma+imaginary)-ndtr((-.5-mu)/sigma+imaginary)
        factor=np.exp(-2j*np.pi*k*mu-2*np.pi**2*sigma**2*k**2)*mass
        fourier+=sign*np.prod(factor,axis=-1)
    norm2=0.
    for mu,sa in zip(centers,signs):
        for nu,sb in zip(centers,signs):
            midpoint=(mu+nu)/2;s=sigma/np.sqrt(2)
            mass=np.prod(ndtr((.5-midpoint)/s)-ndtr((-.5-midpoint)/s))
            overlap=np.exp(-np.sum((mu-nu)**2)/(4*sigma*sigma))/(4*np.pi*sigma*sigma)**1.5
            norm2+=sa*sb*overlap*mass
    if not np.isfinite(fourier).all():raise ValueError('Complex Gaussian integral exceeded numerical range')
    return fourier,float(norm2)


def fourier_field_norm2(k,coefficients,chunk=128):
    """Integral of Re(sum c_j exp(2 pi i k_j x)) squared over the unit cube."""
    k=np.asarray(k).reshape(-1,3);c=np.asarray(coefficients).ravel();total=0.;absolute=0.
    for start in range(0,len(k),chunk):
        ki=k[start:start+chunk];ci=c[start:start+chunk]
        minus=np.ones((len(ki),len(k)));plus=minus.copy()
        for axis in range(3):
            minus*=np.sinc(ki[:,axis,None]-k[None,:,axis])
            plus*=np.sinc(ki[:,axis,None]+k[None,:,axis])
        total+=.5*np.real(ci@(minus@c.conj()+plus@c))
        absolute+=.5*abs(ci)@((abs(minus)+abs(plus))@abs(c))
    return float(total),float(absolute)


def continuous_residual_norm(k,ctf,weights,noise_std,centers,signs,width):
    """Norm of the continuous adjoint residual with a roundoff diagnostic pad."""
    n,nq=np.asarray(ctf).shape;w=np.asarray(weights).reshape(n,2*nq)
    coefficients=(w[:,:nq]+1j*w[:,nq:])*ctf/noise_std
    target_fourier,target_norm2=gaussian_target_integrals(k,centers,signs,width)
    cross=float(np.real(np.sum(coefficients*target_fourier.conj())))
    field_norm2,absolute=fourier_field_norm2(k,coefficients)
    residual2=target_norm2-2*cross+field_norm2
    # This protects against ordinary cancellation at these tested scales, but
    # is not a validated interval-arithmetic bound for all floating-point input.
    pad=20*np.finfo(float).eps*len(coefficients.ravel())*(target_norm2+2*abs(cross)+absolute)
    if residual2 < -pad:raise ValueError('Negative residual norm beyond numerical pad')
    return {'residual_norm':float(np.sqrt(max(0,residual2)+pad)),'target_norm':float(np.sqrt(target_norm2)),
            'target_norm2':target_norm2,'field_norm2':field_norm2,'cross_term':cross,
            'residual_norm2_unpadded':residual2,'squared_norm_roundoff_pad':float(pad)}


def cell_target_coefficients(box,centers,signs,width):
    """Exact coefficients in the orthonormal piecewise-constant cube basis."""
    edges=np.linspace(-.5,.5,box+1);result=np.zeros((box,)*3)
    for mu,sign in zip(np.asarray(centers),signs):
        axes=[np.diff(ndtr((edges-mu[a])/width)) for a in range(3)]
        result+=sign*axes[2][:,None,None]*axes[1][None,:,None]*axes[0][None,None,:]
    return box**1.5*result.ravel()


def cell_forward(k,ctf,coefficients,box,noise_std):
    """Exact Fourier observations of normalized constant cells, not point samples."""
    op=VoxelObservationOperator(k,ctf,box,noise_std*box**1.5);n,nq=op.n,op.q
    values=op.forward(coefficients).reshape(n,2*nq);complex_values=values[:,:nq]+1j*values[:,nq:]
    factor=np.prod(np.sinc(k/box),axis=-1)*np.exp(-1j*np.pi*np.sum(k,axis=-1)/box)
    exact=complex_values*factor
    return np.concatenate([exact.real,exact.imag],axis=1).ravel()


def cell_adjoint(k,ctf,weights,box,noise_std):
    """Exact cell-basis projection of the continuous Fourier adjoint."""
    n,nq=np.asarray(ctf).shape;w=np.asarray(weights).reshape(n,2*nq)
    phase=np.exp(1j*np.pi*np.sum(k,axis=-1)/box)
    c=(w[:,:nq]+1j*w[:,nq:])*phase
    transformed=np.concatenate([c.real,c.imag],axis=1).ravel()
    op=VoxelObservationOperator(k,ctf*np.prod(np.sinc(k/box),axis=-1),box,noise_std*box**1.5)
    return op.adjoint(transformed)


class CellObservationOperator:
    """Cached forward/adjoint for the orthonormal constant-cell cube basis."""
    def __init__(self,k,ctf,box,noise_std):
        self.op=VoxelObservationOperator(k,ctf*np.prod(np.sinc(k/box),axis=-1),box,noise_std*box**1.5)
        self.phase=np.exp(-1j*np.pi*np.sum(k,axis=-1)/box);self.shape=self.op.shape

    def forward(self,coefficients):
        nq=self.op.q;v=self.op.forward(coefficients).reshape(self.op.n,2*nq)
        c=(v[:,:nq]+1j*v[:,nq:])*self.phase
        return np.concatenate([c.real,c.imag],axis=1).ravel()

    def adjoint(self,weights):
        nq=self.op.q;w=np.asarray(weights).reshape(self.op.n,2*nq)
        c=(w[:,:nq]+1j*w[:,nq:])*self.phase.conj()
        return self.op.adjoint(np.concatenate([c.real,c.imag],axis=1).ravel())


class ContinuousObservationGram:
    """Exact continuous observation Gram on the unit cube, at a bounded size.

    Central symmetry makes cosine/sine cross terms zero. Two dense matrices
    store the remaining blocks; no density grid or Fourier bandwidth truncation
    is used to define the unknown class. Cost is quadratic in Fourier samples.
    """
    def __init__(self,k,ctf,noise_std,chunk=128,max_frequencies=8192):
        self.k=np.asarray(k);self.ctf=np.asarray(ctf);self.noise_std=noise_std
        self.n,self.nq=self.ctf.shape;points=self.k.reshape(-1,3);size=len(points)
        if size>max_frequencies:raise ValueError('Dense continuous Gram exceeds explicit frequency budget')
        transfer=self.ctf.ravel()/noise_std
        self.real=np.empty((size,size));self.imag=np.empty((size,size));self.shape=(2*size,2*size)
        for start in range(0,size,chunk):
            ki=points[start:start+chunk];minus=np.ones((len(ki),size));plus=minus.copy()
            for axis in range(3):
                minus*=np.sinc(ki[:,axis,None]-points[None,:,axis]);plus*=np.sinc(ki[:,axis,None]+points[None,:,axis])
            scale=.5*transfer[start:start+chunk,None]*transfer[None]
            self.real[start:start+len(ki)]=scale*(minus+plus);self.imag[start:start+len(ki)]=scale*(minus-plus)
        self.diagonal=self.pack(np.diag(self.real),np.diag(self.imag))

    def pack(self,real,imag):
        return np.concatenate([np.asarray(real).reshape(self.n,self.nq),np.asarray(imag).reshape(self.n,self.nq)],axis=1).ravel()

    def matvec(self,weights):
        w=np.asarray(weights).reshape(self.n,2*self.nq)
        return self.pack(self.real@w[:,:self.nq].ravel(),self.imag@w[:,self.nq:].ravel())

    def target(self,centers,signs,width):
        fourier,norm2=gaussian_target_integrals(self.k,centers,signs,width)
        a=fourier*self.ctf/self.noise_std
        return self.pack(a.real,a.imag),norm2


def continuous_certificate(gram,centers,signs,width,density_radius,alpha=.05,maxiter=100,rtol=.005,cg_rtol=1e-9,callback=None):
    """Classical bias-aware majorization with a feasible continuous dual bound.

    h=ell-A*w belongs to L2, and v=t*h is dual feasible for
    t=min(B/||h||, z/||A h||). Its target pairing is known analytically.
    This certifies the sum-width objective, not optimality of folded-normal q.
    """
    B=float(density_radius);z=norm.ppf(1-alpha/2);a,lnorm2=gram.target(centers,signs,width)
    if B<=0 or lnorm2<=0 or not 0<alpha<1:raise ValueError('Positive radius, nonzero target and valid alpha required')
    base=B*np.sqrt(lnorm2);eps=base*1e-7;w=np.zeros_like(a);gw=np.zeros_like(a);lower=0.;history=[]
    best={'weights':w.copy(),'objective':float(base),'bias':float(base),'noise_sd':0.,'half_width':float(base),'target_norm':float(np.sqrt(lnorm2))}
    error_bound=getattr(gram,'quadrature_error',lambda weights:{'squared_field_norm':0.,'gram_action_norm':0.})
    for it in range(maxiter):
        errors=error_bound(w)
        sd=np.linalg.norm(w);h2=max(0.,lnorm2-2*w@a+w@gw+errors['squared_field_norm']);hnorm=np.sqrt(h2)
        ridge=(z*z/np.hypot(z*sd,eps))/(B*B/np.hypot(B*hnorm,eps))
        system=LinearOperator(gram.shape,matvec=lambda x:gram.matvec(x)+ridge*x,dtype=float)
        pre=gram.preconditioner(ridge) if hasattr(gram,'preconditioner') else LinearOperator(gram.shape,matvec=lambda x:x/(gram.diagonal+ridge),dtype=float)
        cg_count=[0]
        def count_iteration(x):cg_count[0]+=1
        w,info=cg(system,a,x0=w,M=pre,rtol=cg_rtol,atol=0,maxiter=1000,callback=count_iteration);gw=gram.matvec(w)
        errors=error_bound(w);sd=np.linalg.norm(w);h2=lnorm2-2*w@a+w@gw+errors['squared_field_norm']
        pad=50*np.finfo(float).eps*len(w)*(lnorm2+2*abs(w@a)+abs(w@gw))
        if h2 < -pad:raise ValueError('Negative continuous residual norm')
        hnorm=np.sqrt(max(0,h2)+pad);bias=B*hnorm;objective=z*sd+bias
        dual_scale=min(B/max(hnorm,1e-300),z/max(np.linalg.norm(a-gw)+errors['gram_action_norm'],1e-300))
        lower=max(lower,0.,float(dual_scale*(lnorm2-w@a)))
        if objective<best['objective']:
            best={'weights':w.copy(),'objective':float(objective),'bias':float(bias),'noise_sd':float(sd),
                  'half_width':bias_aware_half_width(sd,bias,alpha),'target_norm':float(np.sqrt(lnorm2))}
        gap=max(0,best['objective']-lower)/best['objective']
        history.append({'iteration':it+1,'ridge':float(ridge),'objective':float(objective),'best_objective':best['objective'],
                        'dual_lower_bound':lower,'relative_gap':gap,'cg_info':int(info),'squared_norm_roundoff_pad':float(pad),
                        'quadrature_errors':{k:float(v) for k,v in errors.items()},'cg_iterations':cg_count[0]})
        if callback is not None:callback(history[-1])
        if gap<=rtol:break
    best.update(dual_lower_bound=lower,converged=gap<=rtol,history=history)
    return best
