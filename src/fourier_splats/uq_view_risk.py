"""Classical CVaR and variance comparators for independent view-group samples."""
import numpy as np
from fourier_splats.uq_view_variance import empirical_bernstein_radius


def empirical_tail_bound(values,kappa,epsilon=0.):
    """Minimum of t+kappa*(empirical stop loss + epsilon*(1-t))."""
    x=np.sort(np.asarray(values,float));n=len(x)
    if x.ndim!=1 or n<1 or not np.isfinite(x).all() or np.min(x)<0 or np.max(x)>1 or kappa<1 or epsilon<0:
        raise ValueError('Probability observations, kappa >= 1 and nonnegative CDF error required')
    t=np.unique(np.r_[0.,x,1.]);index=np.searchsorted(x,t,side='right')
    prefix=np.r_[0.,np.cumsum(x)];excess=(prefix[-1]-prefix[index]-(n-index)*t)/n
    objective=t+kappa*(excess+epsilon*(1-t));i=int(np.argmin(objective))
    return float(np.clip(objective[i],0,1)),float(t[i])


def risk_baselines(grouped,center,kappas,delta=.001):
    x=np.asarray(grouped,float)
    if x.ndim!=2 or x.shape[1]!=2 or len(x)<4:raise ValueError('Independent views with two noise groups required')
    mean=x.mean(axis=1);m=len(mean);mu=float(mean.mean())
    rad=empirical_bernstein_radius(mean,0,1,delta/2);lo=max(0,mu-rad);hi=min(1,mu+rad)
    extremum=np.clip(.5,lo,hi);mean_only_variance=float(extremum*(1-extremum))
    joint_rad=empirical_bernstein_radius(mean,0,1,delta/4);joint_lo=max(0,mu-joint_rad);joint_hi=min(1,mu+joint_rad)
    second=(mean-center)**2
    second_upper=float(second.mean())+empirical_bernstein_radius(second,0,max(center**2,(1-center)**2),delta/2)
    unpaired_variance=float(np.clip(second_upper-max(0,joint_lo-center,center-joint_hi)**2,0,.25))
    first,second_half=mean[:m//2],mean[m//2:];epsilon=np.sqrt(np.log(2/delta)/(2*m));rows=[]
    for kappa in kappas:
        dkw,t=empirical_tail_bound(mean,kappa,epsilon)
        split_t=float(np.quantile(first,1-1/kappa,method='inverted_cdf'))
        excess=np.maximum(0,second_half-split_t)
        split_upper=float(excess.mean())+empirical_bernstein_radius(excess,0,1-split_t,delta)
        rows.append(dict(kappa=float(kappa),cvar_dkw=dkw,cvar_dkw_threshold=t,
            cvar_split=min(1.,split_t+kappa*split_upper),cvar_split_threshold=split_t,
            mean_only=min(1.,kappa*hi,hi+np.sqrt((kappa-1)*mean_only_variance)),
            unpaired_variance=min(1.,kappa*joint_hi,joint_hi+np.sqrt((kappa-1)*unpaired_variance))))
    return dict(mean_interval=[lo,hi],mean_only_variance=mean_only_variance,
        unpaired_variance_upper=unpaired_variance,dkw_epsilon=float(epsilon),bounds=rows)


def paired_variance_decomposition(grouped,center,delta=.001):
    x=np.asarray(grouped,float);m=len(x);mean=x.mean(axis=1);mu=float(mean.mean())
    radius=empirical_bernstein_radius(mean,0,1,delta/4);lo=max(0,mu-radius);hi=min(1,mu+radius)
    y=(x[:,0]-center)*(x[:,1]-center);width=max(center**2,(1-center)**2)+center*(1-center)
    log=np.log(4/delta);sampling=float(np.sqrt(2*np.var(y,ddof=1)*log/m));ranging=float(7*width*log/(3*(m-1)))
    subtract=max(0,lo-center,center-hi)**2;raw=float(y.mean())
    return dict(product_estimate=raw,square_root_term=sampling,range_term=ranging,
        distance_subtraction=subtract,unclipped_upper=raw+sampling+ranging-subtract,
        variance_upper=float(np.clip(raw+sampling+ranging-subtract,0,.25)),
        sample_variance_of_group_mean=float(np.var(mean,ddof=1)))
