"""Known bounded-density viewing distributions for simulation controls only."""
import numpy as np
from scipy.spatial.transform import Rotation


def sample_cap_views(count, axis, kappa, rng):
    """Haar conditioned on |R[2,axis]| >= 1-1/kappa.

    A Haar row is uniform on the sphere. Its specified coordinate is uniform
    on [-1,1], so the retained event has probability 1/kappa. The resulting
    density relative to Haar is exactly kappa on that event and zero outside.
    Unused accepted proposals are discarded without further selection.
    """
    if count < 1 or int(count) != count or axis not in (0,1,2) or not np.isfinite(kappa) or kappa < 1:
        raise ValueError('Positive integer count, coordinate axis and kappa >= 1 required')
    pieces=[];remaining=int(count);proposed=accepted=0;cut=1-1/kappa
    while remaining:
        n=max(16,int(np.ceil(remaining*kappa*1.05)))
        rotations=Rotation.random(n,random_state=rng).as_matrix()
        keep=rotations[np.abs(rotations[:,2,axis])>=cut]
        proposed+=n;accepted+=len(keep)
        selected=keep[:remaining];pieces.append(selected);remaining-=len(selected)
    return np.concatenate(pieces),dict(proposed=proposed,accepted=accepted,used=int(count),cut=cut)
