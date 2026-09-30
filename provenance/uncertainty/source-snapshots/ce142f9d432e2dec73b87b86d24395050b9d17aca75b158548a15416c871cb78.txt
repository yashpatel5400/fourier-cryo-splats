"""Per-frequency Taylor disks for continuous Gaussian-orbit likelihood bounds.

Classical Fenchel duality; motivated by focused Fable audit mixture-audit-01.
Complex Fourier errors use two-dimensional disks, not independent rectangles.
All guarantees concern real arithmetic, not validated floating-point intervals.
"""
import numpy as np
from .uq_mixture_refinement import SelectiveCurvatureOrbit


def validate_independent_plane(plane):
    """Reject DC, duplicated frequencies and conjugate pairs in the noise model.

    This checks coordinates only. It cannot establish independence of actual
    image noise, external whitening or CTF estimation from evaluation data.
    """
    p = np.asarray(plane, float)
    if p.ndim != 2 or p.shape[1] != 3 or not len(p) or not np.isfinite(p).all():
        raise ValueError('Nonempty finite three-dimensional frequencies required')
    tolerance = 64*np.finfo(float).eps*max(float(np.max(np.abs(p))), 1.)
    if np.any(np.linalg.norm(p, axis=1) <= tolerance):
        raise ValueError('DC cannot be an independent complex noise coordinate')
    for i in range(len(p)):
        if (np.any(np.linalg.norm(p[i+1:]-p[i], axis=1) <= tolerance)
                or np.any(np.linalg.norm(p[i+1:]+p[i], axis=1) <= tolerance)):
            raise ValueError('Duplicate or conjugate complex noise coordinates')


def disk_residual_dual(residual, jacobian, radii, half_width, *, sweeps=24):
    """Bound min_{|v|<=h, ||e_q||<=eps_q} ||r-Jv-e||^2.

    Real coordinates are packed as all real then all imaginary components.
    Projected gradient merely proposes v; block soft thresholding proposes a
    dual vector. Every saved dual is feasible regardless of convergence.
    """
    r, j, eps, h = (np.asarray(x, float) for x in
                    (residual, jacobian, radii, half_width))
    if (r.ndim != 2 or r.shape[1] % 2 or j.shape != (*r.shape, 3)
            or eps.shape != (len(r), r.shape[1]//2) or h.shape != (3,)
            or np.any(eps < 0) or np.any(h < 0) or sweeps < 0
            or any(not np.isfinite(x).all() for x in (r,j,eps,h))):
        raise ValueError('Compatible finite residual, Jacobian, disks and box required')
    q = r.shape[1]//2
    gram = np.einsum('ndj,ndk->njk', j, j)
    cross = np.einsum('ndj,nd->nj', j, r)
    lipschitz = np.maximum(np.linalg.eigvalsh(gram)[:,-1], 0.)
    v = np.zeros((len(r),3))
    best = np.zeros(len(r)); best_u = np.zeros_like(r)
    primal = np.full(len(r), np.inf)

    def evaluate():
        nonlocal best, best_u, primal
        residual_v = r-np.einsum('ndj,nj->nd',j,v)
        norm = np.hypot(residual_v[:,:q],residual_v[:,q:])
        scale = np.maximum(1.-np.divide(eps,norm,out=np.ones_like(eps),where=norm>0),0.)
        u = residual_v*np.concatenate([scale,scale],axis=1)
        jt = np.einsum('ndj,nd->nj',j,u)
        value = (2*np.sum(r*u,axis=1)-np.sum(u*u,axis=1)
                 -2*np.sum(eps*np.hypot(u[:,:q],u[:,q:]),axis=1)
                 -2*np.sum(h*np.abs(jt),axis=1))
        take = value > best
        best_u[take] = u[take]; best = np.maximum(best,value)
        primal = np.minimum(primal,np.sum(u*u,axis=1))
        return jt

    evaluate()
    # A cheap ordinary box least-squares anchor, followed by disk-aware steps.
    for _ in range(8):
        for axis in range(3):
            gradient = cross[:,axis]-np.sum(gram[:,axis,:]*v,axis=1)
            step = np.divide(gradient,gram[:,axis,axis],out=np.zeros_like(gradient),where=gram[:,axis,axis]>0)
            v[:,axis] = np.clip(v[:,axis]+step,-h[axis],h[axis])
        evaluate()
    for _ in range(sweeps):
        direction = evaluate()
        v = np.clip(v+np.divide(direction,lipschitz[:,None],
                   out=np.zeros_like(direction),where=lipschitz[:,None]>0),-h,h)
    evaluate()
    return dict(lower_squared=best,dual_vectors=best_u,
                primal_squared=primal,gap=np.maximum(primal-best,0.))


class DiskGaussianOrbit(SelectiveCurvatureOrbit):
    """Intersect earlier cell envelopes with the per-frequency disk dual."""
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        validate_independent_plane(self.plane)

    def cell(self,lower,upper):
        original = super().cell(lower,upper)
        half = (np.asarray(upper)-np.asarray(lower))/2
        h = float(half.sum())
        mean,jac,k = self.mean_jacobian(original['middle'])
        displacement = 2*self.kr*np.sin(min(h,np.pi)/2)
        b1,b2 = np.zeros(self.q),np.zeros(self.q)
        for sign in [-1,1]:
            distance = np.linalg.norm(k[:,None,:]+sign*self.centers,axis=-1)
            lo = np.maximum(distance-displacement[:,None],0.)
            hi = distance+displacement[:,None]
            peak = np.clip(self.sigma,lo,hi)
            exponential = np.exp(-.5*(peak/self.sigma)**2)
            b1 += (peak/self.sigma**2*exponential)@self.abs_coefficient
            b2 += ((1+(peak/self.sigma)**2)/self.sigma**2*exponential)@self.abs_coefficient
        radii = .5*h*h*np.abs(self.transfer)*(self.kr**2*b2+self.kr*b1)
        fit = disk_residual_dual(self.observed-mean,jac,radii,half)
        disk_upper = -.5*fit['lower_squared']
        original.update(previous_envelope=original['envelope'].copy(),
            disk_envelope=disk_upper,disk_primal_gap=fit['gap'],
            envelope=np.minimum(original['envelope'],disk_upper))
        if np.any(original['envelope'] < original['exact']-1e-8):
            raise FloatingPointError('Disk envelope below center')
        return original
