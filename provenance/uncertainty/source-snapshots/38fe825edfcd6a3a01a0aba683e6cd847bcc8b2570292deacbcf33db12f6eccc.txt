"""Classical scaled-anchor bounds (Lindsay 1983); no new inference principle."""
import numpy as np
from scipy.special import logsumexp


def scaled_anchor_upper(log_envelopes, log_anchor):
    """Upper bound on max_w sum_i log(sum_c U_ic w_c), for any anchor>0.

    Unlike a feasible mixture, the anchor can be arbitrary. This formula is
    invariant to a common scaling of all anchors. Callers must separately
    establish that U is an upper envelope for every allowed latent value.
    """
    u,z=np.asarray(log_envelopes,float),np.asarray(log_anchor,float)
    if (u.ndim!=2 or min(u.shape)==0 or z.shape!=(u.shape[0],)
            or np.any(np.isnan(u)) or not np.isfinite(z).all()
            or np.any(np.all(np.isneginf(u),axis=1))):
        raise ValueError('Nonempty envelopes and finite compatible anchors required')
    n=len(z)
    return float(z.sum()+n*(np.max(logsumexp(u-z[:,None],axis=0))-np.log(n)))
