"""Matrix-free continuous Fourier Gram with an explicit quadrature remainder.

Tensor Gauss--Legendre quadrature is applied only to finite Fourier sums, not
to arbitrary unknown density. The unknown class remains the full cube L2 ball.
The stated remainder is a real-arithmetic quadrature bound; FINUFFT roundoff is
checked numerically and is not an interval-arithmetic guarantee.
"""
import sys
import numpy as np
import finufft
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator
from .uq_continuous import gaussian_target_integrals


class QuadratureObservationGram:
    def __init__(self,k,ctf,noise_std,order=40,eps=1e-12,preconditioner_rank=256):
        self.k=np.asarray(k);self.ctf=np.asarray(ctf);self.noise_std=noise_std;self.order=order;self.eps=eps
        self.n,self.nq=self.ctf.shape;self.shape=(2*self.n*self.nq,)*2
        self.transfer=self.ctf/noise_std
        nodes,weights=leggauss(order);nodes/=2;weights/=2
        z,y,x=np.meshgrid(nodes,nodes,nodes,indexing='ij')
        self.nodes=[np.ascontiguousarray(a.ravel()) for a in [x,y,z]]
        self.weights=(weights[:,None,None]*weights[None,:,None]*weights[None,None,:]).ravel()
        points=2*np.pi*self.k.reshape(-1,3);self.frequencies=[np.ascontiguousarray(points[:,a]) for a in range(3)]
        self.nthreads=1 if sys.platform=='darwin' else 4
        twice_kernel=np.prod(np.sinc(2*self.k),axis=-1)
        self.diagonal=self.pack(.5*self.transfer**2*(1+twice_kernel),.5*self.transfer**2*(1-twice_kernel))
        maximum=np.max(abs(self.k),axis=(0,1));constant=4*gammaln(order+1)-np.log(2*order+1)-3*gammaln(2*order+1)
        errors=[np.exp(constant+2*order*np.log(4*np.pi*v)) if v else 0. for v in maximum]
        self.kernel_error=float(np.sum(errors))
        self.preconditioner_factors=[];self.preconditioner_diagnostics=[]
        if preconditioner_rank:
            points=self.k.reshape(-1,3);transfer=self.transfer.ravel()
            for sign in [1,-1]:
                diagonal=.5*transfer**2*(1+sign*twice_kernel.ravel());initial=float(diagonal.max())
                factor=np.zeros((len(points),min(preconditioner_rank,len(points))))
                used=0
                for column in range(factor.shape[1]):
                    pivot=int(np.argmax(diagonal));value=diagonal[pivot]
                    if value<=1e-12*max(initial,1):break
                    minus=np.prod(np.sinc(points-points[pivot]),axis=1);plus=np.prod(np.sinc(points+points[pivot]),axis=1)
                    entry=.5*transfer*transfer[pivot]*(minus+sign*plus)
                    if column:entry-=factor[:,:column]@factor[pivot,:column]
                    factor[:,column]=entry/np.sqrt(value);diagonal=np.maximum(0,diagonal-factor[:,column]**2);used=column+1
                factor=factor[:,:used]
                eigenvalues,eigenvectors=eigh(factor.T@factor)
                self.preconditioner_factors.append((factor@eigenvectors,np.maximum(eigenvalues,0)))
                self.preconditioner_diagnostics.append({'rank':used,'maximum_residual_diagonal':float(diagonal.max())})

    def pack(self,real,imag):
        return np.concatenate([np.asarray(real).reshape(self.n,self.nq),np.asarray(imag).reshape(self.n,self.nq)],axis=1).ravel()

    def coefficients(self,weights):
        w=np.asarray(weights).reshape(self.n,2*self.nq)
        return (w[:,:self.nq]+1j*w[:,self.nq:])*self.transfer

    def field(self,weights):
        coefficients=np.ascontiguousarray(self.coefficients(weights).ravel())
        return finufft.nufft3d3(*self.frequencies,coefficients,*self.nodes,isign=1,eps=self.eps,nthreads=self.nthreads).real

    def matvec(self,weights):
        field=self.field(weights);weighted=np.ascontiguousarray(field*self.weights,dtype=complex)
        values=finufft.nufft3d3(*self.nodes,weighted,*self.frequencies,isign=-1,eps=self.eps,nthreads=self.nthreads)
        values=values.reshape(self.n,self.nq)*self.transfer
        return self.pack(values.real,values.imag)

    def quadrature_error(self,weights):
        total=float(np.sum(abs(self.coefficients(weights))))
        return {'squared_field_norm':self.kernel_error*total**2,
                'gram_action_norm':np.sqrt(2)*self.kernel_error*np.linalg.norm(self.transfer)*total}

    def target(self,centers,signs,width):
        ft,norm2=gaussian_target_integrals(self.k,centers,signs,width);a=ft*self.transfer
        return self.pack(a.real,a.imag),norm2

    def preconditioner(self,ridge):
        if not self.preconditioner_factors:
            return LinearOperator(self.shape,matvec=lambda x:x/(self.diagonal+ridge),dtype=float)
        def apply(x):
            x=np.asarray(x).reshape(self.n,2*self.nq);out=[]
            for block,(factor,eigenvalues) in zip([x[:,:self.nq].ravel(),x[:,self.nq:].ravel()],self.preconditioner_factors):
                out.append((block-factor@((factor.T@block)/(ridge+eigenvalues)))/ridge)
            return self.pack(*out)
        return LinearOperator(self.shape,matvec=apply,dtype=float)
