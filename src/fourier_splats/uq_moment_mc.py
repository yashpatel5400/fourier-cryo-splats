"""Finite-simulation candidate tests under an explicitly bounded viewing law.

No experimental noise/view calibration is supplied by this module. Polynomial
profiling is floating-point evaluation, not validated interval arithmetic.
"""
import numpy as np
from scipy.stats import beta,binom


def noisy_amplitude_polynomial(means,noise,triads,direction,noise_variance=1.):
    """Ascending cubic coefficients of the debiased contrast f(a*mean+noise)."""
    m=np.atleast_2d(np.asarray(means,complex));e=np.atleast_2d(np.asarray(noise,complex))
    triads=np.asarray(triads,int).reshape(-1,3);w=np.asarray(direction,float)
    n=m.shape[1];t=len(triads)
    if e.shape!=m.shape or w.shape!=(n+2*t,) or not np.isfinite(m).all() or not np.isfinite(e).all():
        raise ValueError('Matching finite means/noise and contrast required')
    p=w[:n]/2;c=np.zeros((len(m),4))
    c[:,0]=np.abs(e)**2@p-noise_variance*w[:n].sum()
    c[:,1]=2*(m*e.conj()).real@p;c[:,2]=np.abs(m)**2@p
    if t:
        a,b,k=triads.T;weight=(w[n:n+t]-1j*w[n+t:])/2
        ma,mb,mc=m[:,a],m[:,b],m[:,k].conj()
        ea,eb,ec=e[:,a],e[:,b],e[:,k].conj()
        for j,value in enumerate([ea*eb*ec,
            ma*eb*ec+ea*mb*ec+ea*eb*mc,
            ma*mb*ec+ma*eb*mc+ea*mb*mc,ma*mb*mc]):
            c[:,j]+=(value@weight).real
    return c


def cubic_interval_maximum(coefficients,lower,upper):
    c=np.asarray(coefficients,float)
    if c.ndim!=2 or c.shape[1]!=4 or not np.isfinite(c).all() or not np.isfinite([lower,upper]).all() or lower>upper:
        raise ValueError('Finite cubics and ordered finite endpoints required')
    def evaluate(x):return c[:,0]+x*(c[:,1]+x*(c[:,2]+x*c[:,3]))
    arg=np.full(len(c),float(lower));best=evaluate(arg)
    def consider(root):
        valid=(root>=lower)&(root<=upper);value=evaluate(np.where(valid,root,lower))
        take=valid&(value>best);best[take]=value[take];arg[take]=root[take]
    consider(np.full(len(c),float(upper)))
    a,b,d=3*c[:,3],2*c[:,2],c[:,1]
    discriminant=b*b-4*a*d;quadratic=(a!=0)&(discriminant>=0)
    # Stable quadratic formula, retaining the linear and double-zero cases.
    z=np.zeros(len(c));z[quadratic]=-.5*(b[quadratic]+np.copysign(np.sqrt(discriminant[quadratic]),b[quadratic]))
    r=np.full(len(c),np.nan);r[quadratic]=z[quadratic]/a[quadratic];consider(r)
    r=np.full(len(c),np.nan);good=quadratic&(z!=0);r[good]=d[good]/z[good];consider(r)
    r=np.full(len(c),np.nan);linear=(a==0)&(b!=0);r[linear]=-d[linear]/b[linear];consider(r)
    return best,arg


def binomial_upper(count,total,delta):
    if not (isinstance(count,(int,np.integer)) and isinstance(total,(int,np.integer))
            and 0<=count<=total and total>0 and 0<delta<1):
        raise ValueError('Valid integer binomial count and tail probability required')
    return 1. if count==total else float(beta.ppf(1-delta,count+1,total-count))


def binomial_interval(count,total,alpha=.05):
    return (0. if count==0 else float(beta.ppf(alpha/2,count,total-count+1)),
            1. if count==total else float(beta.ppf(1-alpha/2,count+1,total-count)))


def rejection_count(total,probability,level):
    """Smallest integer k with Pr[Bin(total,probability)>=k]<=level."""
    if not (total>0 and 0<=probability<=1 and 0<level<1):raise ValueError('Invalid binomial test')
    lo,hi=0,int(total)+1
    while lo<hi:
        middle=(lo+hi)//2
        if binom.sf(middle-1,total,probability)<=level:hi=middle
        else:lo=middle+1
    return lo


def projected_rejection_probability(count,total,particles,critical):
    """Pointwise 95% interval from independent single-particle simulations.

    These are binomial-law projections, not observed replicate dataset power.
    """
    lower,upper=binomial_interval(count,total)
    return dict(single_particle_probability=count/total,probability_interval=[lower,upper],
        rejection_probability=float(binom.sf(critical-1,particles,count/total)),
        rejection_interval=[float(binom.sf(critical-1,particles,lower)),float(binom.sf(critical-1,particles,upper))])
