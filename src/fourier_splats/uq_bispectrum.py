"""Translation-invariant moment diagnostics; no continuous-pose certificate.

Noise coordinates have independent N(0, variance) real and imaginary parts.
Triads use three distinct nonredundant complex coordinates a,b,c with qa+qb=qc.
"""
import time
import numpy as np
from scipy.sparse import coo_matrix


def frequency_triads(q, maximum=512, seed=261010):
    q = np.asarray(q, float)
    if q.ndim != 2 or q.shape[1] != 2 or not np.isfinite(q).all():
        raise ValueError('Finite two-dimensional frequencies required')
    lookup = {tuple(x): i for i, x in enumerate(q)}
    if len(lookup) != len(q) or (0., 0.) in lookup or any(tuple(-x) in lookup for x in q):
        raise ValueError('Distinct, nonzero, nonopposite frequencies required')
    triples = []
    for a in range(len(q)):
        for b in range(a + 1, len(q)):
            c = lookup.get(tuple(q[a] + q[b]))
            if c is not None and c not in (a, b):
                triples.append((a, b, c))
    triples = np.asarray(triples, int).reshape(-1, 3)
    if maximum is not None and len(triples) > maximum:
        selected = np.random.default_rng(seed).choice(len(triples), maximum, replace=False)
        triples = triples[np.sort(selected)]
    return triples


def moment_features(values, triads, noise_variance=0.):
    """Power and real/imaginary bispectrum, each divided by two.

The divisor makes each coordinate have unit variance at pure unit Gaussian
noise, but does NOT whiten their covariance at nonzero signal.
"""
    values = np.asarray(values, complex)
    triads = np.asarray(triads, int).reshape(-1, 3)
    power = (abs(values) ** 2 - 2 * noise_variance) / 2
    if not len(triads):
        return power
    a, b, c = triads.T
    bispectrum = values[..., a] * values[..., b] * values[..., c].conj() / 2
    return np.concatenate([power, bispectrum.real, bispectrum.imag], axis=-1)


class MomentContrast:
    """Exact Gaussian mean/variance for a fixed cubic contrast.

Uses the finite Gaussian Hermite expansion; includes all cross-feature
covariances. Sparse Hessian assembly avoids a dense third derivative tensor.
"""
    def __init__(self, frequencies, triads, coefficients):
        self.n = int(frequencies)
        self.triads = np.asarray(triads, int).reshape(-1, 3)
        w = np.asarray(coefficients, float)
        n, t = self.n, len(self.triads)
        if w.shape != (n + 2*t,) or not np.isfinite(w).all():
            raise ValueError('One finite weight per moment coordinate required')
        if t and (np.min(self.triads) < 0 or np.max(self.triads) >= n
                  or any(len(set(row)) != 3 or row[0] >= row[1] for row in self.triads)
                  or len(set(map(tuple, self.triads))) != t):
            raise ValueError('Distinct valid triads required')
        self.coefficients = w
        self.power = w[:n] / 2
        self.cubic = (w[n:n+t] - 1j*w[n+t:]) / 2
        rows, cols, values, lookup = [], [], [], {}
        # Each unique off-diagonal real Hessian entry is a linear function
        # of the mean. Entries sharing a pair are added before squaring.
        for triad, weight in zip(self.triads, self.cubic):
            for left, right, third in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
                for il in [0, 1]:
                    for ir in [0, 1]:
                        pair = tuple(sorted((triad[left]+il*n, triad[right]+ir*n)))
                        row = lookup.setdefault(pair, len(lookup))
                        fl = (1 if il == 0 else (-1j if left == 2 else 1j))
                        fr = (1 if ir == 0 else (-1j if right == 2 else 1j))
                        for it in [0, 1]:
                            ft = (1 if it == 0 else (-1j if third == 2 else 1j))
                            rows.append(row); cols.append(triad[third]+it*n)
                            values.append(float(np.real(weight*fl*fr*ft)))
        self.hessian_map = coo_matrix((values, (rows, cols)),
                                     shape=(len(lookup), 2*n)).tocsr()

    def variance_polynomial(self, means, noise_variance=1., batch=128):
        """Coefficients of Var[contrast(a*mean+noise)] in ascending powers of a."""
        means = np.atleast_2d(np.asarray(means, complex))
        v = float(noise_variance)
        if means.shape[1] != self.n or not np.isfinite(means).all() or not np.isfinite(v) or v < 0:
            raise ValueError('Matching finite means and nonnegative variance required')
        coefficients = np.zeros((len(means), 5))
        n = self.n
        for begin in range(0, len(means), batch):
            m = means[begin:begin+batch]
            real = np.concatenate([m.real, m.imag], axis=1)
            power_gradient = real * np.tile(2*self.power, 2)
            gradient = np.zeros_like(real)
            if len(self.triads):
                a, b, c = self.triads.T
                for index, factor, conjugated in [
                    (a, self.cubic*m[:, b]*m[:, c].conj(), False),
                    (b, self.cubic*m[:, a]*m[:, c].conj(), False),
                    (c, self.cubic*m[:, a]*m[:, b], True)]:
                    for row in range(len(m)):
                        np.add.at(gradient[row], index, factor[row].real)
                        np.add.at(gradient[row], index+n,
                                  factor[row].imag * (1 if conjugated else -1))
            off_diagonal = self.hessian_map @ real.T
            rows = coefficients[begin:begin+len(m)]
            rows[:,0] = v*v*4*np.sum(self.power**2)+v**3*4*np.sum(abs(self.cubic)**2)
            rows[:,2] = v*np.sum(power_gradient**2,axis=1)+v*v*np.sum(off_diagonal**2,axis=0)
            rows[:,3] = 2*v*np.sum(power_gradient*gradient,axis=1)
            rows[:,4] = v*np.sum(gradient**2,axis=1)
        return coefficients

    def mean_variance(self, means, noise_variance=1., batch=128):
        means = np.atleast_2d(np.asarray(means, complex))
        expectation = moment_features(means, self.triads) @ self.coefficients
        return expectation, self.variance_polynomial(means,noise_variance,batch).sum(axis=1)


def amplitude_maximum(coefficients, lower, upper):
    """Maximize c0+c2*a²+c3*a³+c4*a⁴ on a closed nonnegative interval.

Check both endpoints and every real stationary point. Floating-point
evaluation of an analytic formula, not outward-rounded interval arithmetic.
"""
    c = np.asarray(coefficients,float)
    if c.ndim != 2 or c.shape[1] != 5 or not np.isfinite(c).all() or np.any(c[:,1] != 0) or not 0 <= lower <= upper:
        raise ValueError('Finite quartics without linear terms and a nonnegative interval required')
    def evaluate(a): return c[:,0]+a*a*(c[:,2]+a*(c[:,3]+a*c[:,4]))
    left=np.full(len(c),float(lower));right=np.full(len(c),float(upper))
    best=evaluate(left);arg=left.copy();value=evaluate(right)
    take=value>best;best[take]=value[take];arg[take]=right[take]
    a,b,d=4*c[:,4],3*c[:,3],2*c[:,2]
    discriminant=b*b-4*a*d
    candidates=[]
    nonzero=a!=0;real=nonzero&(discriminant>=0)
    for sign in [-1.,1.]:
        root=np.full(len(c),np.nan)
        root[real]=(-b[real]+sign*np.sqrt(discriminant[real]))/(2*a[real])
        candidates.append(root)
    root=np.full(len(c),np.nan);linear=(~nonzero)&(b!=0)
    root[linear]=-d[linear]/b[linear];candidates.append(root)
    for root in candidates:
        valid=(root>=lower)&(root<=upper)
        value=evaluate(np.where(valid,root,lower));take=valid&(value>best)
        best[take]=value[take];arg[take]=root[take]
    return best,arg


def project_simplex(x):
    x = np.asarray(x, float)
    u = np.sort(x)[::-1]
    css = np.cumsum(u)-1
    good = u-css/np.arange(1, len(u)+1) > 0
    k = np.flatnonzero(good)[-1]
    z = np.maximum(x-css[k]/(k+1), 0.)
    return z/z.sum()


def moment_hull_projection(atoms, target, max_iterations=2500, maximum_seconds=90.,
                           distance_tolerance=1e-5):
    """Feasible mixture and separating-hyperplane distance bracket.

Accelerated projected gradient with checked backtracking; stopping bounds
do not assume optimizer convergence. All bounds are floating-point diagnostics.
"""
    a, target = np.asarray(atoms, float), np.asarray(target, float)
    if a.ndim != 2 or target.shape != a.shape[1:] or not np.isfinite(a).all() or not np.isfinite(target).all():
        raise ValueError('Finite atom matrix and matching target required')
    start = time.perf_counter()
    w = np.full(len(a), 1/len(a)); y = w.copy(); acceleration = 1.
    # Power iteration estimates a step size; backtracking, not this estimate,
    # verifies each accepted quadratic step.
    z = np.ones(a.shape[1]); z /= np.linalg.norm(z)
    for _ in range(20):
        z = a.T@(a@z)
        length = np.linalg.norm(z)
        if length == 0: break
        z /= length
    lipschitz = max(float(np.sum((a@z)**2)), 1e-12)
    history = []; converged = False; status = 'iteration_limit'; backtracks = 0
    best_w = w.copy(); best_upper = np.inf; best_lower = 0.; best_direction = np.zeros_like(target)
    iteration = 0
    for iteration in range(max_iterations+1):
        if iteration % 25 == 0 or iteration == max_iterations:
            approximation = w@a; residual = target-approximation
            upper = float(np.linalg.norm(residual))
            direction = residual/upper if upper > 0 else np.zeros_like(residual)
            lower = max(0., float(direction@target-np.max(a@direction)))
            if upper < best_upper: best_upper = upper; best_w = w.copy()
            if lower > best_lower: best_lower = lower; best_direction = direction.copy()
            if best_lower > best_upper+1e-8*max(1., best_upper):
                raise ArithmeticError('Hull weak duality failed')
            history.append(dict(iteration=iteration,upper=best_upper,lower=best_lower,
                                gap=best_upper-best_lower,seconds=time.perf_counter()-start))
            if best_upper-best_lower <= distance_tolerance:
                converged = True; status = 'distance_gap'; break
            if time.perf_counter()-start > maximum_seconds:
                status = 'time_limit'; break
        if iteration == max_iterations: break
        residual = y@a-target; gradient = a@residual; fy = .5*(residual@residual)
        for _ in range(60):
            proposal = project_simplex(y-gradient/lipschitz)
            diff = proposal-y; new_residual = proposal@a-target
            fn = .5*(new_residual@new_residual)
            upper_model = fy+gradient@diff+.5*lipschitz*(diff@diff)
            if fn <= upper_model+1e-12*max(1., abs(fn), abs(fy)): break
            lipschitz *= 2; backtracks += 1
        else: raise ArithmeticError('Backtracking did not find a valid step')
        next_acceleration = (1+np.sqrt(1+4*acceleration**2))/2
        y = proposal+(acceleration-1)/next_acceleration*(proposal-w)
        w = proposal; acceleration = next_acceleration
    return dict(weights=best_w,approximation=best_w@a,direction=best_direction,
                distance_upper=best_upper,distance_lower=best_lower,
                gap=max(0.,best_upper-best_lower),converged=converged,status=status,
                iterations=iteration,backtracks=backtracks,history=history,
                seconds=time.perf_counter()-start)
