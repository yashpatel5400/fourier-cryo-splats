"""Skip identically unused NUFFT channels when every rotation radius is zero.

The field still has the original twenty-column lift, spatial quadrature and
random-probe dimension. Translation-only coefficients have spatial degree zero.
"""
import numpy as np
import finufft
from .uq_pose_operator import PolynomialPoseFieldOperator


class ShiftAwarePolynomialPoseFieldOperator(PolynomialPoseFieldOperator):
    def _constant_channel(self, strengths, source, destination):
        if self.backend == 'nufft':
            return finufft.nufft3d3(*source, np.ascontiguousarray(strengths[0]), *destination,
                                    isign=1, eps=self.eps, nthreads=self.nthreads)
        output = np.empty(len(destination[0]), complex)
        s = np.stack(source, axis=1); d = np.stack(destination, axis=1)
        for start in range(0, len(d), 1024):
            output[start:start+1024] = strengths[0]@np.exp(1j*s@d[start:start+1024].T)
        return output

    def _transform_forward(self, strengths):
        if np.any(self.angle != 0):
            return super()._transform_forward(strengths)
        if np.any(strengths[1:] != 0):
            raise AssertionError('Zero-angle operator acquired a nonconstant spatial coefficient')
        output = np.zeros((10, len(self.xyz)), complex)
        output[0] = self._constant_channel(strengths, self.source, self.destination)
        return output

    def _transform_adjoint(self, strengths):
        if np.any(self.angle != 0):
            return super()._transform_adjoint(strengths)
        # Nonconstant moments are discarded only because every coefficient
        # multiplying them in rmatvec/weight_gradient is identically zero.
        output = np.zeros((10, len(self.frequency_points)), complex)
        output[0] = self._constant_channel(strengths, self.destination, self.source)
        return output
