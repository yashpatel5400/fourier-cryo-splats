"""Second-order log-kernel enclosures for the continuous orbit prototype."""
from itertools import product
import numpy as np
from .uq_continuous_mixture import FourierGaussianOrbit


def euler_second_derivatives(angles):
    """Second derivatives of the same Rz(a) Ry(b) Rz(c) row-vector action."""
    a,b,c=np.asarray(angles,float)
    def z(t):
        co,si=np.cos(t),np.sin(t)
        return [np.array([[co,-si,0.],[si,co,0.],[0.,0.,1.]]),
                np.array([[-si,-co,0.],[co,-si,0.],[0.,0.,0.]]),
                np.array([[-co,si,0.],[-si,-co,0.],[0.,0.,0.]])]
    co,si=np.cos(b),np.sin(b)
    y=[np.array([[co,0.,si],[0.,1.,0.],[-si,0.,co]]),
       np.array([[-si,0.,co],[0.,0.,0.],[-co,0.,-si]]),
       np.array([[-co,0.,-si],[0.,0.,0.],[si,0.,-co]])]
    factors=[z(a),y,z(c)]
    answer=np.zeros((3,3,3,3))
    for j,k in product(range(3),repeat=2):
        order=[int(i==j)+int(i==k) for i in range(3)]
        answer[j,k]=factors[0][order[0]]@factors[1][order[1]]@factors[2][order[2]]
    return answer


def quadratic_box_upper(gradient,hessian,half_width,sweeps=12):
    """Concave tangent upper, valid for any anchor in real arithmetic."""
    g,h=np.asarray(gradient,float),np.asarray(hessian,float)
    box=np.asarray(half_width,float)
    if g.shape[-1]!=3 or h.shape!=(*g.shape,3) or box.shape!=(3,) or np.any(box<0):
        raise ValueError('Compatible three-dimensional quadratic box required')
    h=.5*(h+np.swapaxes(h,-1,-2))
    eigen=np.linalg.eigvalsh(h)
    # A numerical guard, not an interval-arithmetic certificate.
    pad=32*np.finfo(float).eps*np.maximum(np.max(np.abs(eigen),axis=-1),1.)
    shift=np.maximum(eigen[...,-1],0.)+pad
    q=shift[...,None,None]*np.eye(3)-h
    v=np.zeros_like(g)
    def bound():
        qv=np.einsum('...jk,...k->...j',q,v)
        return (.5*np.sum(v*qv,axis=-1)+np.sum(box*np.abs(g-qv),axis=-1)
                +.5*shift*np.sum(box**2))
    upper=bound()
    for _ in range(sweeps):
        for j in range(3):
            grad=g[...,j]-np.sum(q[...,j,:]*v,axis=-1)
            step=np.divide(grad,q[...,j,j],out=np.zeros_like(grad),where=q[...,j,j]>0)
            v[...,j]=np.clip(v[...,j]+step,-box[j],box[j])
        upper=np.minimum(upper,bound())
    return upper


class CurvatureGaussianOrbit(FourierGaussianOrbit):
    def log_kernel_derivatives(self,angles):
        from .uq_continuous_mixture import euler_rotation_jacobian
        mean,jac,k=self.mean_jacobian(angles)
        _,dr=euler_rotation_jacobian(angles)
        d2r=euler_second_derivatives(angles)
        dk=np.einsum('qi,aij->aqj',self.plane,dr)
        d2k=np.einsum('qi,abij->abqj',self.plane,d2r)
        size=len(self.centers)
        coefficient=self.coefficients[:size]+1j*self.coefficients[size:]
        gradient=np.zeros((self.q,3),complex);hessian=np.zeros((self.q,3,3),complex)
        for sign,c in [(-1,coefficient),(1,coefficient.conj())]:
            d=k[:,None,:]+sign*self.centers
            weight=np.exp(-.5*np.sum(d*d,axis=-1)/self.sigma**2)*c[None,:]
            gradient-=np.einsum('qmi,qm->qi',d,weight/self.sigma**2)
            hessian+=np.einsum('qmi,qmj,qm->qij',d,d,weight/self.sigma**4)
            hessian-=np.sum(weight/self.sigma**2,axis=1)[:,None,None]*np.eye(3)
        second=np.einsum('aqi,qij,bqj->qab',dk,hessian,dk)+np.einsum('qi,abqi->qab',gradient,d2k)
        second=self.transfer[:,:,None,None]*second[None,:,:,:]
        second=np.concatenate([second.real,second.imag],axis=1)
        residual=self.observed-mean
        log_gradient=np.einsum('nd,ndj->nj',residual,jac)
        log_hessian=np.einsum('nd,ndjk->njk',residual,second)-np.einsum('ndj,ndk->njk',jac,jac)
        return -.5*np.sum(residual**2,axis=1),log_gradient,log_hessian,k,residual

    def cell(self,lower,upper):
        original=super().cell(lower,upper)
        half=(np.asarray(upper)-np.asarray(lower))/2;h=float(half.sum())
        log0,g,hess,k,residual=self.log_kernel_derivatives(original['middle'])
        displacement=2*self.kr*np.sin(min(h,np.pi)/2)
        b1,b2,b3=(np.zeros(self.q) for _ in range(3))
        for sign in [-1,1]:
            distance=np.linalg.norm(k[:,None,:]+sign*self.centers,axis=-1)
            lo=np.maximum(distance-displacement[:,None],0.);hi=distance+displacement[:,None]
            peak=np.clip(self.sigma,lo,hi)
            exponential=np.exp(-.5*(peak/self.sigma)**2)
            b1+=(peak/self.sigma**2*exponential)@self.abs_coefficient
            b2+=((1+(peak/self.sigma)**2)/self.sigma**2*exponential)@self.abs_coefficient
            peak=np.clip(3**.25*self.sigma,lo,hi)
            b3+=((peak**3/self.sigma**6+3*peak/self.sigma**4)
                 *np.exp(-.5*(peak/self.sigma)**2))@self.abs_coefficient
        transfer=np.abs(self.transfer)
        m1=transfer*(self.kr*b1)
        m2=transfer*(self.kr**2*b2+self.kr*b1)
        m3=transfer*(self.kr**3*b3+3*self.kr**2*b2+self.kr*b1)
        complex_residual=residual[:,:self.q]+1j*residual[:,self.q:]
        third=h**3/6*np.sum((np.abs(complex_residual)+h*m1)*m3+3*m1*m2,axis=1)
        quadratic=quadratic_box_upper(g,hess,half)
        improved=np.minimum(log0+quadratic+third,0.)
        original.update(first_order_envelope=original['envelope'].copy(),
            quadratic_envelope=improved,third_log_remainder=third,
            quadratic_increment_upper=quadratic,
            envelope=np.minimum(original['envelope'],improved))
        if np.any(original['envelope']<log0-1e-8):
            raise FloatingPointError('Curvature envelope below center')
        return original
