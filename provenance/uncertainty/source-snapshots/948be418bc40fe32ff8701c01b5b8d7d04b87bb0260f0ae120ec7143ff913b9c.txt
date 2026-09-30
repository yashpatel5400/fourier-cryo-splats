"""Continuous SO(3) likelihood envelopes for Hermitian Fourier Gaussians.

Known noise/CTFs, zero translation, one fixed homogeneous map. The statistical
dual is classical (Lindsay, 1983); these are floating-point implementations of
real-arithmetic enclosures, not validated interval arithmetic.
"""
from itertools import product
import time
import numpy as np
from scipy.special import logsumexp
from .uq_physics import evaluate_pairs
from .uq_mixture_validation import mixture_envelope_upper


def euler_rotation_jacobian(angles):
    """Rz(a) Ry(b) Rz(c) and its three derivatives (radians)."""
    a, b, c = np.asarray(angles, float)
    def z(t):
        co, si = np.cos(t), np.sin(t)
        return (np.array([[co,-si,0.],[si,co,0.],[0.,0.,1.]]),
                np.array([[-si,-co,0.],[co,-si,0.],[0.,0.,0.]]))
    co, si = np.cos(b), np.sin(b)
    ry = np.array([[co,0.,si],[0.,1.,0.],[-si,0.,co]])
    dy = np.array([[-si,0.,co],[0.,0.,0.],[-co,0.,-si]])
    ra, da = z(a); rc, dc = z(c)
    return ra@ry@rc, np.stack([da@ry@rc, ra@dy@rc, ra@ry@dc])


def linear_box_residual_dual(residual_squared, cross, gram, half_width,
                             sweeps=8):
    """Lower squared residual for min_|v|<=h ||r-Jv||², via feasible duals.

    cross=J^T r, gram=J^T J. Coordinate descent merely supplies dual anchors.
    No convergence or primal optimality is needed for the lower bound.
    """
    r2, c, g = map(lambda x: np.asarray(x, float),
                    (residual_squared, cross, gram))
    h = np.asarray(half_width, float)
    if (c.shape[-1] != len(h) or g.shape != (*c.shape, len(h))
            or r2.shape != c.shape[:-1] or np.any(h < 0)):
        raise ValueError('Incompatible nonnegative box and residual arrays')
    v = np.zeros_like(c)
    def dual():
        gv = np.einsum('...jk,...k->...j', g, v)
        return r2 - np.sum(v*gv, axis=-1) - 2*np.sum(h*np.abs(c-gv), axis=-1)
    lower = np.maximum(dual(), 0.)
    for _ in range(sweeps):
        for j in range(len(h)):
            derivative = c[..., j] - np.sum(g[..., j, :]*v, axis=-1)
            step = np.divide(derivative, g[..., j, j],
                             out=np.zeros_like(derivative), where=g[..., j, j] > 0)
            v[..., j] = np.clip(v[..., j]+step, -h[j], h[j])
        lower = np.maximum(lower, dual())
    return np.maximum(lower, 0.)


class FourierGaussianOrbit:
    """A fixed map observed through known image-specific transfer functions.

    All coordinates/coefficients are as in uq_physics.evaluate_pairs. Transfer
    already includes supplied noise whitening; observed real/imag noise has
    unit variance. No DC or conjugate duplicate should be supplied as an
    independent noisy coordinate by the caller.
    """
    def __init__(self, plane, centers, sigma, coefficients, transfer, observed):
        self.plane = np.asarray(plane, float)
        self.centers = np.asarray(centers, float)
        self.sigma = np.broadcast_to(np.asarray(sigma, float), (len(self.centers),))
        self.coefficients = np.asarray(coefficients, float)
        self.transfer = np.asarray(transfer, float)
        self.observed = np.asarray(observed, float)
        q = len(self.plane)
        if (self.plane.shape != (q,3) or self.centers.ndim != 2
                or self.centers.shape[1] != 3 or np.any(self.sigma <= 0)
                or self.coefficients.shape != (2*len(self.centers),)
                or self.transfer.ndim != 2 or self.transfer.shape[1] != q
                or self.observed.shape != (len(self.transfer),2*q)
                or any(not np.isfinite(x).all() for x in [self.plane,
                    self.centers,self.sigma,self.coefficients,self.transfer,self.observed])):
            raise ValueError('Finite compatible Gaussian orbit inputs required')
        self.n, self.q = self.transfer.shape
        self.normalization = -self.n*self.q*np.log(2*np.pi)
        self.kr = np.linalg.norm(self.plane, axis=1)
        self.center_norm = np.linalg.norm(self.centers,axis=1)
        self.abs_coefficient = np.abs(self.coefficients[:len(self.centers)]
                                     +1j*self.coefficients[len(self.centers):])
        self.calls = 0

    def mean_jacobian(self, angles):
        rotation, derivative = euler_rotation_jacobian(angles)
        k = self.plane@rotation
        value, gradient = evaluate_pairs(k,self.centers,self.sigma,
                                          self.coefficients,gradient=True)
        dk = np.einsum('qi,aij->aqj', self.plane, derivative)
        # dk[a,q,j] = sum_i plane[q,i] dR[a,i,j].
        df = np.einsum('qj,aqj->qa', gradient, dk)
        mean = self.transfer*value[None,:]
        jac = self.transfer[:,:,None]*df[None,:,:]
        return (np.concatenate([mean.real,mean.imag],axis=1),
                np.concatenate([jac.real,jac.imag],axis=1), k)

    def log_kernel(self, angles):
        mean, _, _ = self.mean_jacobian(angles)
        return -.5*np.sum((self.observed-mean)**2,axis=1)

    def cell(self, lower, upper):
        lower, upper = np.asarray(lower,float), np.asarray(upper,float)
        if (lower.shape != (3,) or upper.shape != (3,) or np.any(upper < lower)
                or not np.isfinite(lower).all() or not np.isfinite(upper).all()):
            raise ValueError('Finite ordered three-dimensional Euler box required')
        middle, half = (lower+upper)/2, (upper-lower)/2
        mean, jac, k = self.mean_jacobian(middle)
        residual = self.observed-mean
        r2 = np.sum(residual**2,axis=1)
        h = float(half.sum())
        displacement = 2*self.kr*np.sin(min(h,np.pi)/2)
        b1, b2 = np.zeros(self.q), np.zeros(self.q)
        for sign in [-1,1]:
            distance = np.linalg.norm(k[:,None,:]+sign*self.centers,axis=-1)
            lo = np.maximum(distance-displacement[:,None],0.)
            hi = distance+displacement[:,None]
            peak = np.clip(self.sigma,lo,hi)
            exponential = np.exp(-.5*(peak/self.sigma)**2)
            b1 += (peak/self.sigma**2*exponential)@self.abs_coefficient
            b2 += ((1+(peak/self.sigma)**2)/self.sigma**2*exponential)@self.abs_coefficient
        abs_transfer = np.abs(self.transfer)
        remainder = .5*h*h*np.linalg.norm(abs_transfer*(self.kr**2*b2+self.kr*b1),axis=1)
        mean_radius = h*np.linalg.norm(abs_transfer*(self.kr*b1),axis=1)
        cross = np.einsum('ndj,nd->nj',jac,residual)
        gram = np.einsum('ndj,ndk->njk',jac,jac)
        linear_lower = linear_box_residual_dual(r2,cross,gram,half)
        residual_lower = np.maximum(np.sqrt(linear_lower)-remainder,0.)
        residual_lower = np.maximum(residual_lower,np.sqrt(r2)-mean_radius)
        # Independent per-frequency rectangular enclosure. On a rotation
        # orbit, |k| is fixed: each Hermitian pair is an even real sum and an
        # odd imaginary difference of two exponentials in dot(k, center).
        # This also keeps high-frequency noise coordinates from being spent
        # freely inside one large image-wide mean ball.
        maximum_dot=self.kr[:,None]*self.center_norm[None,:]
        cosine=np.divide(k@self.centers.T,maximum_dot,
                         out=np.zeros_like(maximum_dot),where=maximum_dot>0)
        angle=np.arccos(np.clip(cosine,-1.,1.));rho=min(h,np.pi)
        dot_lo=maximum_dot*np.cos(np.minimum(angle+rho,np.pi))
        dot_hi=maximum_dot*np.cos(np.maximum(angle-rho,0.))
        radial=self.kr[:,None]**2+self.center_norm[None,:]**2
        def pair(dot):
            # max(0,...) removes harmless cancellation in squared distances.
            plus=np.exp(-.5*np.maximum(radial-2*dot,0.)/self.sigma**2)
            minus=np.exp(-.5*np.maximum(radial+2*dot,0.)/self.sigma**2)
            return plus+minus,plus-minus
        re_lo=pair(np.clip(np.zeros_like(dot_lo),dot_lo,dot_hi))[0]
        re_hi=pair(np.maximum(np.abs(dot_lo),np.abs(dot_hi)))[0]
        im_lo=pair(dot_lo)[1];im_hi=pair(dot_hi)[1]
        c=self.coefficients;size=len(self.centers)
        def sum_interval(lo,hi,coefficient):
            return (lo@np.maximum(coefficient,0.)+hi@np.minimum(coefficient,0.),
                    hi@np.maximum(coefficient,0.)+lo@np.minimum(coefficient,0.))
        re_lo,re_hi=sum_interval(re_lo,re_hi,c[:size])
        im_lo,im_hi=sum_interval(im_lo,im_hi,c[size:])
        transfer=np.concatenate([self.transfer,self.transfer],axis=1)
        bound_a=transfer*np.concatenate([re_lo,im_lo])[None,:]
        bound_b=transfer*np.concatenate([re_hi,im_hi])[None,:]
        coordinate_distance=np.maximum(np.maximum(np.minimum(bound_a,bound_b)-self.observed,
                            self.observed-np.maximum(bound_a,bound_b)),0.)
        coordinate_lower=np.linalg.norm(coordinate_distance,axis=1)
        residual_lower=np.maximum(residual_lower,coordinate_lower)
        exact, envelope = -.5*r2, -.5*residual_lower**2
        if np.any(envelope < exact-1e-8):
            raise FloatingPointError('Cell upper below its center likelihood')
        self.calls += 1
        return dict(lower=lower,upper=upper,middle=middle,
                    exact=exact,envelope=envelope,remainder=remainder,
                    coordinate_residual_lower=coordinate_lower)


def continuous_mixture_upper(orbit, *, initial_bins=(4,2,4), max_splits=1024,
                             refit_every=64, mixture_iterations=500,
                             tolerance=1., wall_seconds=600., callback=None):
    """Bracket continuous likelihood; retain cover and best replayable anchors.

    The returned likelihoods include Gaussian normalization. Every complete
    covering partition yields a valid real-arithmetic upper bound, including
    the initial coarse cover if no splits/iterations are allowed.
    """
    if (len(initial_bins)!=3 or any(int(b)!=b or b<1 for b in initial_bins)
            or max_splits<0 or refit_every<1 or mixture_iterations<0
            or tolerance<=0 or wall_seconds<=0):
        raise ValueError('Positive budgets/bins and nonnegative split limit required')
    begin=time.perf_counter()
    endpoints=[np.linspace(0.,end,int(b)+1) for b,end in zip(initial_bins,[2*np.pi,np.pi,2*np.pi])]
    leaves={};support=[];support_angles=[];history=[];split_history=[]
    def add_cell(low,high):
        cell=orbit.cell(low,high);index=len(support)
        support.append(cell['exact']);support_angles.append(cell['middle'])
        leaves[index]=cell
        return index
    for ind in product(*(range(int(b)) for b in initial_bins)):
        add_cell([endpoints[j][i] for j,i in enumerate(ind)],
                 [endpoints[j][i+1] for j,i in enumerate(ind)])
    best_lower,best_upper=-np.inf,np.inf
    best_anchor=None;best_weights=None;best_support_count=0;best_cover=None
    best_fit_weights=None;best_fit_count=0
    log_z=None
    status='split_limit'
    for step in range(max_splits+1):
        refitted=step%refit_every==0 or step==max_splits
        if refitted:
            fitted=mixture_envelope_upper(np.array(support).T,
                    tolerance=min(.01,tolerance/10),max_iterations=mixture_iterations)
            log_z=logsumexp(np.array(support).T+np.log(fitted.weights)[None,:],axis=1)
            primal=float(log_z.sum())
            if primal>best_lower:
                best_lower=primal;best_fit_weights=fitted.weights.copy();best_fit_count=len(support)
        ids=list(leaves)
        uppers=np.array([leaves[i]['envelope'] for i in ids])
        log_scores=logsumexp(uppers-log_z[None,:],axis=1)
        worst=int(np.argmax(log_scores))
        upper=float(log_z.sum()+orbit.n*(log_scores[worst]-np.log(orbit.n)))
        if upper<best_lower-1e-7:
            raise FloatingPointError('Continuous mixture bracket crossed')
        if upper<best_upper:
            best_upper=upper;best_anchor=log_z.copy()
            best_weights=fitted.weights.copy();best_support_count=len(fitted.weights)
            best_cover=dict(ids=np.array(ids),envelopes=uppers.copy(),
                lower=np.array([leaves[i]['lower'] for i in ids]),
                upper=np.array([leaves[i]['upper'] for i in ids]))
        gap=max(best_upper-best_lower,0.)
        if refitted or gap<=tolerance or time.perf_counter()-begin>=wall_seconds:
            row=dict(splits=step,leaves=len(leaves),support_count=len(support),
                     lower=best_lower+orbit.normalization,upper=best_upper+orbit.normalization,
                     gap=gap,seconds=time.perf_counter()-begin,
                     finite_fit_gap=fitted.dual_gap,finite_fit_converged=fitted.converged)
            history.append(row)
            if callback is not None:callback(row)
        if gap<=tolerance:
            status='converged';break
        if time.perf_counter()-begin>=wall_seconds:
            status='wall_limit';break
        if step==max_splits:break
        parent_id=ids[worst];parent=leaves[parent_id]
        axis=int(np.argmax(parent['upper']-parent['lower']))
        mid=(parent['lower'][axis]+parent['upper'][axis])/2
        left_upper=parent['upper'].copy();left_upper[axis]=mid
        right_lower=parent['lower'].copy();right_lower[axis]=mid
        left=add_cell(parent['lower'],left_upper)
        right=add_cell(right_lower,parent['upper'])
        del leaves[parent_id]
        split_history.append([parent_id,axis,left,right])
    # The best upper's cover may predate final leaves. Retain its exact value;
    # the final refined cover with the same anchors also upper-bounds the
    # likelihood, but interval formulas need not be numerically monotone.
    final_ids=np.array(list(leaves))
    final_envelopes=np.array([leaves[i]['envelope'] for i in final_ids])
    replay_upper=float(best_anchor.sum()+orbit.n*(
        np.max(logsumexp(final_envelopes-best_anchor[None,:],axis=1))-np.log(orbit.n)))
    if replay_upper<best_upper:
        best_upper=replay_upper
        best_cover=dict(ids=final_ids.copy(),envelopes=final_envelopes.copy(),
            lower=np.array([leaves[i]['lower'] for i in final_ids]),
            upper=np.array([leaves[i]['upper'] for i in final_ids]))
    if best_upper-best_lower<=tolerance:
        status='converged'
    return dict(log_likelihood_lower=best_lower+orbit.normalization,
        log_likelihood_upper=best_upper+orbit.normalization,
        gap=max(best_upper-best_lower,0.),status=status,
        converged=best_upper-best_lower<=tolerance,seconds=time.perf_counter()-begin,
        history=history,split_history=np.asarray(split_history,dtype=int).reshape(-1,4),
        leaf_ids=final_ids,leaf_lower=np.array([leaves[i]['lower'] for i in final_ids]),
        leaf_upper=np.array([leaves[i]['upper'] for i in final_ids]),
        leaf_log_envelopes=final_envelopes,support_angles=np.array(support_angles),
        support_log_kernels=np.array(support),best_anchor_log_z=best_anchor,
        best_cover_ids=best_cover['ids'],best_cover_log_envelopes=best_cover['envelopes'],
        best_cover_lower=best_cover['lower'],best_cover_upper=best_cover['upper'],
        best_anchor_weights=best_weights,best_anchor_support_count=best_support_count,
        best_feasible_weights=best_fit_weights,best_feasible_support_count=best_fit_count,
        replay_upper_on_final_cover=replay_upper+orbit.normalization,
        normalization=orbit.normalization,
        scope='Known unit noise, fixed CTF, zero shifts and fixed Gaussian map; floating-point real-arithmetic bounds.')
