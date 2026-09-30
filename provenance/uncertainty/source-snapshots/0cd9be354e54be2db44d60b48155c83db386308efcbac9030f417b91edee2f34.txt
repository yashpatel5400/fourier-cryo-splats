"""Noise-metric design transforms, with no claim that the metric is calibrated.

The metric can be learned from an independent training pool. Final feature
variance must be audited separately, without selecting weights on that audit
pool. This module only reparameterizes a continuous fixed-pose design problem.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator


class CenteredNoiseMetric:
    """Symmetric powers of O_i C O_i.T, where O_i is a Fourier centering phase."""
    def __init__(self,covariance_proxy,phases):
        self.covariance=np.asarray(covariance_proxy,float);self.phases=np.asarray(phases,complex)
        if self.phases.ndim!=2 or not np.isfinite(self.phases).all() or not np.allclose(abs(self.phases),1.):
            raise ValueError('Finite unit-modulus phase matrix required')
        self.n,self.nq=self.phases.shape;self.m=2*self.nq
        if self.covariance.shape!=(self.m,self.m) or not np.isfinite(self.covariance).all() or not np.allclose(self.covariance,self.covariance.T):
            raise ValueError('Finite symmetric covariance proxy required')
        self.values,self.vectors=eigh(self.covariance)
        if self.values[0]<=0:raise ValueError('Positive definite covariance proxy required')
        self.inverse_sqrt_norm=float(self.values[0]**-.5)

    def apply(self,vector,power=-.5):
        """Accept one stacked vector or a matrix with stacked vectors as columns."""
        x=np.asarray(vector,float);matrix=x.ndim==2
        if x.ndim not in [1,2] or x.shape[0]!=self.n*self.m:raise ValueError('Wrong stacked observation dimension')
        x=x.reshape(self.n,self.m,-1)
        z=(x[:,:self.nq]+1j*x[:,self.nq:])*self.phases.conjugate()[:,:,None]
        raw=np.concatenate([z.real,z.imag],axis=1)
        transform=(self.vectors*self.values**power)@self.vectors.T
        result=np.einsum('ab,nbk->nak',transform,raw)
        z=(result[:,:self.nq]+1j*result[:,self.nq:])*self.phases[:,:,None]
        out=np.concatenate([z.real,z.imag],axis=1).reshape(self.n*self.m,-1)
        return out if matrix else out[:,0]


class NoiseMetricGram:
    """Transform a continuous Gram by W G W, W=C_centered^(-1/2).

    Its standard Euclidean noise norm is a DESIGN PROXY. The returned interval
    from a generic optimizer is not calibrated for the actual noise law.
    Quadrature error still protects the continuous density residual exactly as
    in the underlying implementation (subject to its numerical caveats).
    """
    def __init__(self,gram,metric):
        self.gram=gram;self.metric=metric;self.shape=gram.shape
        if self.shape[0]!=metric.n*metric.m:raise ValueError('Incompatible metric and Gram')
        factors=[]
        for block,(f,_) in enumerate(gram.preconditioner_factors):
            packed=np.zeros((gram.n,2*gram.nq,f.shape[1]))
            packed[:,block*gram.nq:(block+1)*gram.nq,:]=f.reshape(gram.n,gram.nq,-1)
            factors.append(metric.apply(packed.reshape(self.shape[0],-1)))
        if factors:
            h=np.concatenate(factors,axis=1);values,vectors=eigh(h.T@h)
            self.factor=h@vectors;self.values=np.maximum(values,0.)
        else:self.factor=None

    def target(self,centers,signs,width):
        a,norm2=self.gram.target(centers,signs,width)
        return self.metric.apply(a),norm2

    def matvec(self,x):return self.metric.apply(self.gram.matvec(self.metric.apply(x)))

    def quadrature_error(self,x):
        error=self.gram.quadrature_error(self.metric.apply(x))
        return {'squared_field_norm':error['squared_field_norm'],
                'gram_action_norm':self.metric.inverse_sqrt_norm*error['gram_action_norm']}

    def preconditioner(self,ridge):
        if ridge<=0:raise ValueError('Positive ridge required')
        if self.factor is None:return LinearOperator(self.shape,matvec=lambda x:x/ridge,dtype=float)
        h=self.factor;values=self.values
        return LinearOperator(self.shape,matvec=lambda x:(x-h@((h.T@x)/(ridge+values)))/ridge,dtype=float)


def shrunk_second_moment(observations,shrinkage=.2):
    """Positive empirical design metric; arbitrary signal means are retained."""
    y=np.asarray(observations,float)
    if y.ndim!=2 or min(y.shape)<1 or not np.isfinite(y).all() or not 0<shrinkage<=1:
        raise ValueError('Finite nonempty training matrix and shrinkage in (0,1] required')
    moment=y.T@y/len(y);scale=float(np.trace(moment)/y.shape[1])
    if scale<=0:raise ValueError('Positive empirical energy required')
    return (1-shrinkage)*moment+shrinkage*scale*np.eye(y.shape[1])
