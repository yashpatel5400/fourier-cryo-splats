"""Noise-replicated nuisance envelopes and a viewing-variance bound.

All guarantees concern a specified simulator. No experimental covariance,
viewing-density or amplitude calibration is inferred by these routines.
"""
import numpy as np


def cubic_cell_events(coefficients,lower,upper,cells,threshold):
    """Whether each cubic exceeds threshold anywhere in each amplitude cell."""
    c=np.asarray(coefficients,float)
    if (c.ndim!=2 or c.shape[1]!=4 or not np.isfinite(c).all()
        or not lower<upper or cells<1 or not np.isfinite(threshold)):
        raise ValueError('Finite cubics, threshold and nonempty cell partition required')
    edges=np.linspace(lower,upper,cells+1)
    value=c[:,0,None]+edges*(c[:,1,None]+edges*(c[:,2,None]+edges*c[:,3,None]))
    events=np.maximum(value[:,:-1],value[:,1:])>threshold
    a,b,d=3*c[:,3],2*c[:,2],c[:,1];discriminant=b*b-4*a*d
    good=(a!=0)&(discriminant>=0);z=np.zeros(len(c))
    z[good]=-.5*(b[good]+np.copysign(np.sqrt(discriminant[good]),b[good]))
    roots=[];root=np.full(len(c),np.nan);root[good]=z[good]/a[good];roots.append(root)
    root=np.full(len(c),np.nan);valid=good&(z!=0);root[valid]=d[valid]/z[valid];roots.append(root)
    root=np.full(len(c),np.nan);valid=(a==0)&(b!=0);root[valid]=-d[valid]/b[valid];roots.append(root)
    for root in roots:
        index=np.flatnonzero((root>=lower)&(root<=upper));x=root[index]
        score=c[index,0]+x*(c[index,1]+x*(c[index,2]+x*c[index,3]))
        hit=score>threshold;index=index[hit];x=x[hit]
        cell=np.minimum(cells-1,np.floor((x-lower)*cells/(upper-lower)).astype(int))
        events[index,cell]=True
    return events


def grouped_amplitude_envelope(coefficients,lower,upper,cells,threshold):
    """Return two conditional-noise group envelopes per independent view.

    Shape: (views, 2, replicas, 4). Grouped mean-then-maximum is bounded
    above by the old per-noise maximum-then-mean. Its conditional expectation
    still dominates every actual amplitude's event probability.
    """
    c=np.asarray(coefficients,float)
    if c.ndim!=4 or c.shape[1]!=2 or c.shape[-1]!=4 or c.shape[2]<1:
        raise ValueError('Two noise groups per view, with cubic coefficients')
    events=cubic_cell_events(c.reshape(-1,4),lower,upper,cells,threshold)
    events=events.reshape(c.shape[:-1]+(cells,))
    grouped=events.mean(axis=2).max(axis=2)
    individual=events.max(axis=3).mean(axis=2)
    return grouped,individual


def empirical_bernstein_radius(values,lower,upper,delta):
    """One-sided Maurer--Pontil (2009), Thm 4, scaled to [lower,upper].

    Calling this with delta/2 on both signs gives a two-sided interval with
    total failure delta. The sample unit must be an independent view, not a
    correlated conditional-noise replica.
    """
    x=np.asarray(values,float)
    if (x.ndim!=1 or len(x)<2 or not np.isfinite(x).all() or not 0<delta<1
        or lower>upper or np.min(x)<lower-1e-14 or np.max(x)>upper+1e-14):
        raise ValueError('Independent bounded observations and valid tail probability required')
    log=np.log(2/delta);variance=float(np.var(x,ddof=1));n=len(x)
    return float(np.sqrt(2*variance*log/n)+7*(upper-lower)*log/(3*(n-1)))


def viewing_probability_bounds(grouped,individual,center,kappas,delta=.001):
    """Three separate, prespecified bounds using the same simulation budget.

    Conditional independence of columns, given view, is essential for the
    product moment. The old independent calibration supplies the fixed center.
    Each method separately has failure at most delta; no free minimum across
    the three methods or simultaneous-error claim is made.
    """
    x=np.asarray(grouped,float);raw=np.asarray(individual,float)
    if x.ndim!=2 or x.shape[1]!=2 or raw.shape!=x.shape or not 0<=center<=1:
        raise ValueError('Two paired noise groups and an independent center required')
    if not np.isfinite(x).all() or np.min(x)<0 or np.max(x)>1 or np.any(x>raw+1e-14) or np.max(raw)>1:
        raise ValueError('Valid probability envelopes required')
    mean=x.mean(axis=1);raw_mean=raw.mean(axis=1)
    ordinary_upper=min(1.,float(mean.mean())+empirical_bernstein_radius(mean,0.,1.,delta))
    individual_upper=min(1.,float(raw_mean.mean())+empirical_bernstein_radius(raw_mean,0.,1.,delta))
    # Allocate delta/2 to two-sided mean confidence and delta/2 to the
    # one-sided centered second moment; dependence between these is allowed.
    radius=empirical_bernstein_radius(mean,0.,1.,delta/4)
    lower=max(0.,float(mean.mean())-radius);upper=min(1.,float(mean.mean())+radius)
    product=(x[:,0]-center)*(x[:,1]-center)
    y_lower=-center*(1-center);y_upper=max(center**2,(1-center)**2)
    product_upper=float(product.mean())+empirical_bernstein_radius(product,y_lower,y_upper,delta/2)
    distance=max(0.,lower-center,center-upper)
    variance_upper=min(.25,max(0.,product_upper-distance**2))
    bounds=[]
    for kappa in kappas:
        if kappa<1 or not np.isfinite(kappa):raise ValueError('Finite viewing ratio at least one required')
        covariance_bound=min(1.,kappa*upper,upper+np.sqrt((kappa-1)*variance_upper))
        bounds.append(dict(kappa=float(kappa),individual_ratio=min(1.,kappa*individual_upper),
            grouped_ratio=min(1.,kappa*ordinary_upper),view_variance=float(covariance_bound)))
    return dict(views=len(x),center=float(center),grouped_mean=float(mean.mean()),
        individual_mean=float(raw_mean.mean()),ordinary_grouped_upper=ordinary_upper,
        ordinary_individual_upper=individual_upper,joint_mean_interval=[lower,upper],
        centered_product_mean=float(product.mean()),centered_product_upper=product_upper,
        view_variance_upper=variance_upper,bounds=bounds)
