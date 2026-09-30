"""Optional reusable CPU NUFFT plans for unchanged Fourier/quadrature grids.

These subclasses preserve the frozen reference implementations. They change
transform setup reuse, not integration nodes, arithmetic precision or bounds.
Grids must remain unchanged after construction; weights may change freely.
Official API: https://finufft.readthedocs.io/en/latest/python.html
"""
import numpy as np
import finufft
from .uq_continuous_quadrature import QuadratureObservationGram
from .uq_pose_operator import PolynomialPoseFieldOperator


class _PlanCache:
    def _execute_type3(self,name,source,destination,strengths,isign):
        strengths=np.ascontiguousarray(strengths,dtype=np.complex128)
        transforms=1 if strengths.ndim==1 else strengths.shape[0]
        parameters=(transforms,isign,self.eps,self.nthreads)
        if not hasattr(self,'_transform_plans'):self._transform_plans={}
        cached=self._transform_plans.get(name)
        if cached is None or cached[0]!=parameters:
            plan=finufft.Plan(3,3,n_trans=transforms,isign=isign,eps=self.eps,
                              dtype='complex128',nthreads=self.nthreads)
            plan.setpts(*source,*destination)
            self._transform_plans[name]=(parameters,plan)
        return self._transform_plans[name][1].execute(strengths)


class PlannedQuadratureObservationGram(_PlanCache,QuadratureObservationGram):
    def field(self,weights):
        return self._execute_type3('field',self.frequencies,self.nodes,
                                  self.coefficients(weights).ravel(),1).real

    def matvec(self,weights):
        field=self.field(weights)
        values=self._execute_type3('adjoint',self.nodes,self.frequencies,field*self.weights,-1)
        values=values.reshape(self.n,self.nq)*self.transfer
        return self.pack(values.real,values.imag)


class PlannedPolynomialPoseFieldOperator(_PlanCache,PolynomialPoseFieldOperator):
    def _transform_forward(self,strengths):
        if self.backend!='nufft':return super()._transform_forward(strengths)
        return self._execute_type3('forward',self.source,self.destination,strengths,1)

    def _transform_adjoint(self,strengths):
        if self.backend!='nufft':return super()._transform_adjoint(strengths)
        return self._execute_type3('adjoint',self.destination,self.source,strengths,1)
