"""Continuous Hilbert-space projection anchored at the target functional.

This projection is an inner optimization relaxation only. It does not redefine
the unknown density class. Exact sinc/Gaussian kernel identities construct the
coordinates in real arithmetic; floating-point orthogonality is checked, not
certified by interval arithmetic.
"""
import numpy as np


def enrich_target_projection(projection,gram,weights,rtol=1e-10):
    """Add the omitted part of A*weights using continuous Gram products.

    Since ell is already in the basis, the new target coordinate is zero.
    The resulting projection represents A*weights exactly in real arithmetic.
    A materially negative residual is an error, never silently square-rooted.
    gram must represent the SAME continuous observation inner product used
    to construct projection. This is a numerical Hilbert Gram-Schmidt step,
    subject to the module's floating-point orthogonality caveat.
    """
    f=np.asarray(projection['factor'],float);b=np.asarray(projection['target'],float)
    w=np.asarray(weights,float)
    if w.shape!=(len(f),) or not np.isfinite(w).all():raise ValueError('Finite compatible weights required')
    gw=gram.matvec(w);coordinates=f.T@w
    residual_column=gw-f@coordinates
    squared=float(w@residual_column);scale=max(float(abs(w@gw)),float(coordinates@coordinates),1.)
    if squared < -1e-8*scale:raise FloatingPointError('Adaptive projection lost PSD substantially')
    diagnostic={'omitted_field_norm2':squared,'relative_omitted_field_norm2':squared/scale,'added':False}
    if squared>rtol*scale:
        projection['factor']=np.column_stack([f,residual_column/np.sqrt(squared)])
        projection['target']=np.append(b,0.)
        diagnostic['added']=True
        projection['diagnostics']['rank']=projection['factor'].shape[1]
    projection['diagnostics'].setdefault('adaptive_enrichments',[]).append(diagnostic)
    return diagnostic


def target_anchored_projection(gram,target_vector,target_norm2,rank=256,rtol=1e-11):
    """Pivoted Hilbert Gram-Schmidt: first ell, then observation kernels.

    Returns F_ij=<A_i*,phi_j> and b_j=<ell,phi_j>, where phi_j span the target
    and selected observation representers. Thus ||b-F.T w|| <= ||ell-A*w||.
    Full-weight-space dual v=sum_j u_j phi_j has A v=F u and norm(v)=norm(u).
    """
    a=np.asarray(target_vector,float);m=gram.shape[0];n,nq=gram.n,gram.nq
    if a.shape!=(m,) or target_norm2<=0 or rank<1 or not np.isfinite(a).all():raise ValueError('Valid target and positive rank required')
    limit=min(int(rank),m+1);f=np.zeros((m,limit));f[:,0]=a/np.sqrt(target_norm2)
    b=np.zeros(limit);b[0]=np.sqrt(target_norm2)
    diagonal=np.asarray(gram.diagonal,float).copy()-f[:,0]**2
    initial=max(float(np.max(gram.diagonal)),1.);minimum=float(diagonal.min())
    if minimum < -1e-8*initial:raise ValueError('Invalid target/observation joint Gram')
    diagonal=np.maximum(diagonal,0.);pivots=[];used=1
    points=gram.k.reshape(-1,3);transfer=gram.transfer.ravel()
    index=np.arange(m).reshape(n,2*nq);real_index=index[:,:nq].ravel();imag_index=index[:,nq:].ravel()
    for j in range(1,limit):
        pivot=int(np.argmax(diagonal));value=diagonal[pivot]
        if value<=rtol*initial:break
        particle,part=divmod(pivot,2*nq);frequency=part%nq;is_real=part<nq;point=particle*nq+frequency
        minus=np.prod(np.sinc(points-points[point]),axis=1);plus=np.prod(np.sinc(points+points[point]),axis=1)
        block=.5*transfer*transfer[point]*(minus+(1 if is_real else -1)*plus)
        column=np.zeros(m);column[real_index if is_real else imag_index]=block
        column-=f[:,:j]@f[pivot,:j];f[:,j]=column/np.sqrt(value)
        residual=diagonal-f[:,j]**2;minimum=min(minimum,float(residual.min()))
        if minimum < -1e-7*initial:raise FloatingPointError('Projection lost PSD substantially')
        diagonal=np.maximum(residual,0.);pivots.append(pivot);used=j+1
    return {'factor':f[:,:used],'target':b[:used],
            'diagnostics':{'rank':used,'requested_rank':rank,'pivots':pivots,
              'maximum_observation_residual_diagonal':float(diagonal.max()),
              'minimum_unclipped_residual_diagonal':minimum,'relative_stop_tolerance':rtol,
              'scope':'Exact-kernel continuous projection in real arithmetic; inner relaxation, not a replacement for the ambient final audit'}}
