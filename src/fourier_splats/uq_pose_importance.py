"""Classical defensive ACG proposal on SO(3), for integration diagnostics."""
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp


class PoseProposal:
    """Antipodally symmetric quaternion mixture; density wrt normalized Haar."""
    def __init__(self, centers, covariances, weights, haar_weight=.05):
        self.centers=Rotation.from_quat(np.asarray(centers,float))
        self.covariances=np.asarray(covariances,float)
        self.weights=np.asarray(weights,float)
        self.haar_weight=float(haar_weight)
        n=len(self.weights)
        if self.covariances.shape!=(n,3,3) or len(self.centers)!=n or np.any(self.weights<=0) or not np.isclose(self.weights.sum(),1) or not 0<haar_weight<1:
            raise ValueError('Positive compatible mixture weights, covariances and defensive mass required')
        self.local_covariance=self.covariances/4
        self.cholesky=np.linalg.cholesky(self.local_covariance)
        self.precision=np.linalg.inv(self.local_covariance)
        self.logdet=np.linalg.slogdet(self.local_covariance)[1]

    def log_density(self, quaternions):
        samples=Rotation.from_quat(np.asarray(quaternions,float))
        parts=[np.full(len(samples),np.log(self.haar_weight))]
        for i in range(len(self.weights)):
            relative=(self.centers[i].inv()*samples).as_quat()
            quadratic=np.einsum('ni,ij,nj->n',relative[:,:3],self.precision[i],relative[:,:3])+relative[:,3]**2
            parts.append(np.log1p(-self.haar_weight)+np.log(self.weights[i])-.5*self.logdet[i]-2*np.log(quadratic))
        return logsumexp(np.array(parts),axis=0)

    def sample(self, count, rng):
        probabilities=np.r_[self.haar_weight,(1-self.haar_weight)*self.weights]
        labels=rng.choice(len(probabilities),size=count,p=probabilities)
        out=np.empty((count,4))
        for label in range(len(probabilities)):
            indices=np.flatnonzero(labels==label)
            raw=rng.normal(size=(len(indices),4))
            if label:
                raw[:,:3]=raw[:,:3]@self.cholesky[label-1].T
            raw/=np.linalg.norm(raw,axis=1)[:,None]
            out[indices]=raw if not label else (self.centers[label-1]*Rotation.from_quat(raw)).as_quat()
        return out,labels


def importance_summary(log_kernels, log_proposal):
    """Direct importance integral estimates and usual diagnostic ESS, no CI."""
    x=np.asarray(log_kernels,float)-np.asarray(log_proposal,float)[:,None]
    total=logsumexp(x,axis=0)
    weights=np.exp(x-total)
    return dict(log_integrals=(total-np.log(len(x))).tolist(),
                log_ratio=float(total[1]-total[0]),
                ess=(1/np.sum(weights**2,axis=0)).tolist(),
                maximum_weight=weights.max(axis=0).tolist())
