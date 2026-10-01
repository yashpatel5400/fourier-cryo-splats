"""FINUFFT plan reuse for a fixed quadrature Gram; same mathematical operator."""
import finufft
import numpy as np


class CachedQuadratureGram:
    """Wrap an existing Gram without mutating its archived implementation.

    Fixed source/target coordinates allow FFT/interpolation planning once per
    geometry. This is FINUFFT's documented plan interface, not a new transform.
    Separate wrappers are needed for changing geometry or concurrent execution.
    """
    def __init__(self, gram, *, eps=None, nthreads=1):
        if nthreads < 1:
            raise ValueError('Positive thread count required')
        self.base = gram
        self.eps = gram.eps if eps is None else eps
        self.nthreads = nthreads
        self._forward_plan = finufft.Plan(3, 3, n_trans=1, eps=self.eps,
            isign=1, dtype='complex128', nthreads=nthreads)
        self._forward_plan.setpts(*gram.frequencies, *gram.nodes)
        self._reverse_plan = finufft.Plan(3, 3, n_trans=1, eps=self.eps,
            isign=-1, dtype='complex128', nthreads=nthreads)
        self._reverse_plan.setpts(*gram.nodes, *gram.frequencies)

    def __getattr__(self, name):
        return getattr(self.base, name)

    def field(self, weights):
        coefficients = np.ascontiguousarray(self.base.coefficients(weights).ravel())
        return self._forward_plan.execute(coefficients).real

    def matvec(self, weights):
        field = self.field(weights)
        values = self._reverse_plan.execute(np.ascontiguousarray(field*self.base.weights, dtype=complex))
        values = values.reshape(self.base.n, self.base.nq)*self.base.transfer
        return self.base.pack(values.real, values.imag)
