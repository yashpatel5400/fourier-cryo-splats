"""Spectral upper bounds retaining the one density error shared by particles.

For columns grouped as H_j and ||v_j|| <= 1, positive d_j give
  ||sum H_j v_j|| <= sqrt(sum d_j) ||[H_j / sqrt(d_j)]||_op.
The constant residual is a singleton block with coefficient exactly one.
Relaxing its coefficient to a unit ball can only increase this bound. The
quadratic pose blocks similarly relax u tensor u to an independent unit ball.
These are classical operator-norm inequalities, not a new general theorem.
"""
import numpy as np
from scipy.linalg import eigvalsh


def grouped_spectral_bound(blocks,weighting='frobenius'):
    """Return a feasible upper bound; eigensolves use ordinary float64 arithmetic."""
    blocks=[np.asarray(a,dtype=np.float64).reshape(len(a),-1) for a in blocks]
    lengths=np.array([np.linalg.norm(a) for a in blocks]);active=lengths>0
    if not active.any():return 0.
    blocks=[a for a,keep in zip(blocks,active) if keep];lengths=lengths[active]
    if weighting=='frobenius':d=lengths
    elif weighting=='uniform':d=np.ones(len(blocks))
    else:raise ValueError('Unknown block weighting')
    matrix=np.concatenate([a/np.sqrt(di) for a,di in zip(blocks,d)],axis=1)
    gram=matrix@matrix.T if matrix.shape[0]<matrix.shape[1] else matrix.T@matrix
    # A small roundoff pad is a numerical precaution, not interval arithmetic.
    top=float(eigvalsh(gram,subset_by_index=[len(gram)-1,len(gram)-1],check_finite=False)[0])
    pad=20*np.finfo(float).eps*max(1,len(gram))*np.linalg.norm(gram,ord=np.inf)
    return float(np.sqrt(d.sum()*max(0,top+pad)))


def shared_density_bounds(residual,linear_blocks,quadratic_blocks):
    """Upper bounds on the density multiplier, before multiplication by B.

quadratic_blocks include their Taylor factor 1/2. They must contain every
ordered derivative pair, so ||u tensor u|| equals ||u|| squared.
"""
    h=np.asarray(residual).reshape(-1,1);linear=list(linear_blocks);quadratic=list(quadratic_blocks)
    triangle=float(np.linalg.norm(h)+sum(np.linalg.norm(a) for a in linear+quadratic))
    first=grouped_spectral_bound(linear);second=grouped_spectral_bound(quadratic)
    separate=float(np.linalg.norm(h)+first+second)
    affine=float(grouped_spectral_bound([h]+linear)+second)
    joint=grouped_spectral_bound([h]+linear+quadratic)
    values={'block_triangle':triangle,'separate_spectral':separate,'affine_linear_spectral':affine,'joint_spectral':joint}
    values['minimum_valid_bound']=min(values.values())
    return values
