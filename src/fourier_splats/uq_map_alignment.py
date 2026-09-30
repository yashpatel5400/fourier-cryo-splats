"""Pilot-only rigid map alignment; a local search, not an accuracy guarantee.

Array coordinates are z,y,x throughout. Returned affine matrices map output
indices to input indices and can be applied unchanged to both native half maps.
"""
import itertools
import numpy as np
from scipy.ndimage import affine_transform
from scipy.optimize import minimize
from scipy.signal import resample
from scipy.spatial.transform import Rotation


def apply_map_alignment(volume, matrix, offset):
    return affine_transform(np.asarray(volume, dtype=float), matrix, offset,
                            order=1, mode='constant', cval=0., prefilter=False)


def _axes(volume):
    weights = np.maximum(volume.ravel(), 0.)
    if weights.sum() <= 0: raise ValueError('Alignment requires positive pilot/source mass')
    weights /= weights.sum()
    xyz = np.indices(volume.shape).reshape(3, -1).T
    center = weights@xyz; delta = xyz-center
    eigenvalues, vectors = np.linalg.eigh((delta.T*weights)@delta)
    return center, vectors[:, ::-1], eigenvalues[::-1]


def _filtered(volume, pixel_size, lowpass_A, box):
    q = np.fft.fftfreq(volume.shape[0], d=pixel_size)
    z, y, x = np.meshgrid(q, q, q, indexing='ij')
    # Smooth radial fourth-order lowpass, followed by Fourier resampling.
    result = np.fft.ifftn(np.fft.fftn(volume)*np.exp(-lowpass_A**4*(x*x+y*y+z*z)**2)).real
    for axis in range(3): result = resample(result, box, axis=axis)
    return result


def align_map_to_pilot(source, pilot, pixel_size, lowpass_A=40., box=32,
                       max_evaluations=300):
    source = np.asarray(source, dtype=float); pilot = np.asarray(pilot, dtype=float)
    if source.shape != pilot.shape or source.ndim != 3 or len(set(source.shape)) != 1:
        raise ValueError('Matched cubic volumes in the same physical field required')
    if not np.isfinite(source).all() or not np.isfinite(pilot).all():
        raise ValueError('Nonfinite map')
    if (not np.isfinite(pixel_size) or not np.isfinite(lowpass_A) or pixel_size <= 0
            or lowpass_A <= 0 or not isinstance(box, int) or not 8 <= box <= source.shape[0]
            or not isinstance(max_evaluations, int) or max_evaluations < 1):
        raise ValueError('Invalid physical or working-grid scale')
    if np.std(source) == 0 or np.std(pilot) == 0: raise ValueError('Constant alignment map')
    moving = _filtered(source, pixel_size, lowpass_A, box)
    fixed = _filtered(pilot, pixel_size, lowpass_A, box)
    sc, sa, se = _axes(moving); tc, ta, te = _axes(fixed)
    fixed0 = fixed-fixed.mean(); fixed_norm = np.linalg.norm(fixed0)
    if fixed_norm <= 0: raise ValueError('Constant pilot')
    starts = []; best = None
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([-1., 1.], repeat=3):
            rotation0 = ta@np.eye(3)[:, perm]@np.diag(signs)@sa.T
            def affine(parameters):
                rotation = Rotation.from_rotvec(parameters[:3]).as_matrix()@rotation0
                matrix = rotation.T
                return matrix, sc-matrix@(tc+parameters[3:])
            def loss(parameters):
                matrix, offset = affine(parameters)
                moved = apply_map_alignment(moving, matrix, offset)
                moved -= moved.mean(); denominator = np.linalg.norm(moved)*fixed_norm
                return 1. if denominator <= 0 else -float(np.sum(moved*fixed0)/denominator)
            initial_loss = loss(np.zeros(6))
            fit = minimize(loss, np.zeros(6), method='Powell',
                bounds=[(-np.pi/3,np.pi/3)]*3+[(-6.,6.)]*3,
                options={'maxfev':max_evaluations,'maxiter':40,'xtol':1e-4,'ftol':1e-6})
            # A bounded local solve can terminate worse than its initial point.
            parameters = fit.x if fit.fun <= initial_loss else np.zeros(6)
            value = min(float(fit.fun), initial_loss); matrix, offset = affine(parameters)
            row = {'start':len(starts),'permutation':list(perm),'signs':list(signs),
                'hand_determinant':int(round(np.linalg.det(matrix))),
                'initial_correlation':-initial_loss,'selected_correlation':-value,
                'local_converged':bool(fit.success),'message':str(fit.message),'evaluations':int(fit.nfev),
                'parameters':parameters.tolist(),'matrix':matrix.tolist(),'working_offset':offset.tolist()}
            starts.append(row)
            if best is None or value < best[0]: best = (value, row)
    row = best[1]; scale = source.shape[0]/box
    matrix = np.array(row['matrix']); offset = np.array(row['working_offset'])*scale
    return {'matrix':matrix,'offset':offset,'selected_start':row['start'], 'starts':starts,
        'source_principal_values':se.tolist(),'pilot_principal_values':te.tolist(),
        'config':{'working_box':box,'pixel_size_A':pixel_size,'lowpass_A':lowpass_A,
                  'max_evaluations_per_start':max_evaluations,'full_map_interpolation_order':1},
        'scope':'Pilot-selected rigid frame and hand from 48 local searches; no reference access or global optimality claim'}
