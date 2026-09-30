"""Weight derivatives for the cubic field; no empirical optimizer is invoked.

The audited cubic operator stays unchanged. These derivatives let a separately
specified design algorithm choose weights before a fresh spectral certificate.
"""
import numpy as np
from .uq_cubic_pose import CubicPoseFieldOperator


class DifferentiableCubicPoseFieldOperator(CubicPoseFieldOperator):
    def moment_weight_gradient(self, column_vector, moments):
        """Gradient in weights of the moment-paired scaled polynomial field."""
        columns = self._to_particle_columns(np.asarray(column_vector)/self.denominators)
        moments = np.asarray(moments)
        if moments.shape != (*self.k.shape[:-1], 20):
            raise ValueError('Twenty moments for each particle frequency required')
        polynomial = np.einsum('nqap,na->nqp', self.polynomials, columns)
        gradient = self.transfer*np.sum(polynomial*moments, axis=-1)
        return np.concatenate([gradient.real, -gradient.imag], axis=1).ravel()

    def weight_gradient(self, column_vector, spatial_vector):
        """Derivative of v.T F(w) u at fixed geometry and column scales."""
        test = np.asarray(spatial_vector, float).reshape(len(self.xyz))
        transformed = self._transform_adjoint((self.monomials*(self.sqrt_quad*test)[None]).astype(complex))
        return self.moment_weight_gradient(column_vector, transformed.T.reshape(self.n, self.nq, 20))
