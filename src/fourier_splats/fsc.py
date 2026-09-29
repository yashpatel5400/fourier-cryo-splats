import numpy as np

def fsc(f1,f2,pixel_size):
    if f1.shape!=f2.shape or f1.ndim!=3: raise ValueError('Matched 3D volumes required')
    d=f1.shape[0]; q=np.arange(-d//2,d//2)
    z,y,x=np.meshgrid(q,q,q,indexing='ij'); rad=np.sqrt(x*x+y*y+z*z)
    shells=np.floor(rad+0.5).astype(int)
    rows=[]
    for s in range(1,d//2-1):
        a=f1[shells==s];b=f2[shells==s]
        den=np.sqrt(np.sum(np.abs(a)**2,dtype=np.float64)*np.sum(np.abs(b)**2,dtype=np.float64))
        corr=np.real(np.sum(a*np.conj(b),dtype=np.complex128))/den if den>0 else np.nan
        rows.append((s,s/(d*pixel_size),corr,len(a)))
    return np.array(rows)

def resolution(curve,threshold=0.143):
    """First crossing sustained for two shells; never extrapolate beyond sampled support."""
    c=np.asarray(curve);r=c[:,2]; freq=c[:,1]
    if not len(c):return {'angstrom':None,'status':'no_shells'}
    if not np.isfinite(r[0]) or r[0]<threshold:return {'angstrom':None,'status':'below_threshold_at_first_shell'}
    for j in range(1,len(c)-1):
        if r[j]<threshold and r[j+1]<threshold:
            prev=j-1
            # Avoid interpolating across any earlier isolated crossing.
            if r[prev]<=threshold: f=freq[j]
            else: f=freq[prev]+(freq[j]-freq[prev])*(r[prev]-threshold)/(r[prev]-r[j])
            return {'angstrom':float(1/f),'status':'crossing','shell':int(c[j,0])}
    return {'angstrom':float(1/freq[-1]),'status':'censored_at_sampled_limit'}
