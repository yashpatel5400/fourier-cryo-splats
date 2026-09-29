"""Feasible nonlinear pose/density stress tests, not global optimality bounds.

For fixed linear-estimator weights and pose u, eliminate the L2 density ball
exactly. Projected gradient ascent then supplies feasible lower bounds on the
largest absolute bias; only an independent certificate supplies an upper bound.
"""
import numpy as np
import torch
from scipy.spatial.transform import Rotation


class TorchPoseAdjoint:
    def __init__(self, operator, q, weights, angle, shift, device='cpu', dtype=torch.float64):
        self.device=device;self.dtype=dtype;self.angle=angle;self.shift=shift
        def tensor(a):return torch.as_tensor(np.asarray(a).copy(),dtype=dtype,device=device)
        op=operator.op;self.k=tensor(op.k);self.q=tensor(q);self.xyz=tensor(operator.xyz)
        self.scale=2*np.pi/op.box;w=np.asarray(weights).reshape(op.n,2*op.q)
        self.wr=tensor(w[:,:op.q]*op.transfer);self.wi=tensor(w[:,op.q:]*op.transfer)

    def __call__(self,u):
        """A(u)'w with right rotation perturbations and detector shifts."""
        v=u[:,:3];x,y,z=v.unbind(-1);zero=torch.zeros_like(x)
        skew=torch.stack([zero,-z,y,z,zero,-x,-y,x,zero],-1).reshape(-1,3,3)
        theta=self.angle*torch.linalg.vector_norm(v,dim=-1)
        rotation=(torch.eye(3,dtype=self.dtype,device=self.device)[None]
                  +self.angle*torch.sinc(theta/torch.pi)[:,None,None]*skew
                  +.5*self.angle**2*torch.sinc(theta/(2*torch.pi))[:,None,None]**2*(skew@skew))
        k=self.k@rotation
        phase=self.scale*(k@self.xyz.T+self.shift*torch.einsum('nqa,na->nq',self.q,u[:,3:])[:,:,None])
        return torch.sum(self.wr[:,:,None]*torch.cos(phase)-self.wi[:,:,None]*torch.sin(phase),dim=(0,1))


def exact_pose_adjoint(operator,q,weights,u,angle,shift):
    """Independent float64 SciPy rotation/direct complex-summation check."""
    op=operator.op;w=np.asarray(weights).reshape(op.n,2*op.q);u=np.asarray(u)
    k=op.k@Rotation.from_rotvec(angle*u[:,:3]).as_matrix()
    phase_shift=shift*np.einsum('nqa,na->nq',q,u[:,3:])
    out=np.zeros(operator.shape[1])
    for i in range(op.n):
        coef=(w[i,:op.q]+1j*w[i,op.q:])*op.transfer[i]
        phase=2j*np.pi/op.box*(k[i]@operator.xyz.T+phase_shift[i,:,None])
        out+=(coef@np.exp(phase)).real
    return out


def eliminated_bias(adjoint,nominal_adjoint,ell,pilot,radius,sign):
    """Worst signed bias over the density ball, for this one feasible pose."""
    h=adjoint-ell
    return sign*torch.dot(pilot,adjoint-nominal_adjoint)+radius*torch.linalg.vector_norm(h)
