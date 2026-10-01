"""Classical continuous quaternion-kernel proposal; no inferential guarantee."""
import numpy as np
from scipy.spatial.transform import Rotation
from .uq_pose_importance import PoseProposal


class CatalogPoseProposal:
    """Haar + local ACG modes + a weighted isotropic ACG rotation catalogue.

    All three parts are normalized relative to normalized Haar measure. Kernel
    centres are fixed before the independent importance samples are generated.
    The catalogue is a proposal, not a discrete replacement for the integral.
    """
    def __init__(self, local, centers, weights, rotation_std, masses=(.05, .45, .5)):
        if not isinstance(local, PoseProposal):
            raise TypeError('Expected a local PoseProposal')
        self.local = local
        self.centers = Rotation.from_quat(np.asarray(centers, float)).as_quat()
        self.weights = np.asarray(weights, float)
        self.masses = np.asarray(masses, float)
        self.s = float(rotation_std)/2
        if (self.weights.shape != (len(self.centers),) or
                not np.all(np.isfinite(self.weights)) or np.any(self.weights < 0) or
                not np.isclose(self.weights.sum(), 1) or self.masses.shape != (3,) or
                not np.all(np.isfinite(self.masses)) or np.any(self.masses <= 0) or
                not np.isclose(self.masses.sum(), 1) or not 0 < self.s <= 1):
            raise ValueError('Normalized weights/masses and a finite kernel scale in (0,1] are required')

    def log_density(self, quaternions, chunk=256):
        samples = Rotation.from_quat(np.asarray(quaternions, float)).as_quat()
        if not isinstance(chunk, int) or chunk < 1:
            raise ValueError('Positive integer chunk required')
        modes_density = np.zeros(len(samples))
        sample_rotations = Rotation.from_quat(samples)
        for i, weight in enumerate(self.local.weights):
            relative = (self.local.centers[i].inv()*sample_rotations).as_quat()
            quadratic = np.einsum('ni,ij,nj->n', relative[:,:3], self.local.precision[i], relative[:,:3])+relative[:,3]**2
            modes_density += weight*np.exp(-.5*self.local.logdet[i])/np.square(quadratic)
        density = self.masses[0] + self.masses[1]*modes_density
        for begin in range(0, len(samples), chunk):
            end = min(begin+chunk, len(samples))
            dot2 = np.square(self.centers @ samples[begin:end].T)
            # Roundoff can place an inner product infinitesimally outside S^3.
            np.clip(dot2, 0, 1, out=dot2)
            quadratic = (1-dot2)/self.s**2 + dot2
            kernel = self.s**-3 / np.square(quadratic)
            density[begin:end] += self.masses[2]*(self.weights @ kernel)
        return np.log(density)

    def sample(self, count, rng):
        if not isinstance(count, int) or count < 1:
            raise ValueError('Positive integer sample count required')
        family = rng.choice(3, size=count, p=self.masses)
        component = np.full(count, -1, dtype=np.int32)
        out = np.empty((count,4))
        for kind in range(3):
            indices = np.flatnonzero(family==kind)
            if not len(indices):
                continue
            raw = rng.normal(size=(len(indices),4))
            if kind==1:
                labels = rng.choice(len(self.local.weights),size=len(indices),p=self.local.weights)
                raw[:,:3] = np.einsum('nij,nj->ni',self.local.cholesky[labels],raw[:,:3])
                centers = self.local.centers[labels]
                component[indices] = labels
            elif kind==2:
                labels = rng.choice(len(self.weights),size=len(indices),p=self.weights)
                raw[:,:3] *= self.s
                centers = Rotation.from_quat(self.centers[labels])
                component[indices] = labels
            raw /= np.linalg.norm(raw,axis=1)[:,None]
            out[indices] = raw if kind==0 else (centers*Rotation.from_quat(raw)).as_quat()
        return out,family,component
