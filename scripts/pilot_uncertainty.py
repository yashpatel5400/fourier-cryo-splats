#!/usr/bin/env python3
"""Development experiment for physical Gaussian features and nuisance bounds.

This deliberately small exact-dictionary experiment checks feasibility, NOT
real-world calibration. The coefficient ball and pose radii are known simulation
assumptions. No inference image or true coefficient tunes an interval.
"""
import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.spatial.transform import Rotation
from scipy.stats import norm

from fourier_splats.uncertainty import optimize_certificate
from fourier_splats.uq_physics import (
    density_functionals, image_design, perturbed_images, pose_jacobian,
    pose_remainder_bounds, pose_design_derivatives, realify,
)

ROOT = Path(__file__).resolve().parents[1]


def geometry(rng, n, preferred=False):
    grid = np.array([(x, y) for y in range(-5, 6) for x in range(-5, 6)
                     if 0 < x*x+y*y <= 25 and (y > 0 or (y == 0 and x > 0))])
    q = np.broadcast_to(grid, (n, len(grid), 2)).copy()
    if preferred:
        # A cone of tilts with free in-plane rotation; orientations are fixed
        # by simulation, so no orientation estimator leaks image noise.
        r = (Rotation.from_rotvec(rng.normal(size=(n, 3))*0.07)
             * Rotation.from_euler('z', rng.uniform(0, 2*np.pi, (n,1)))).as_matrix()
    else:
        r = Rotation.random(n, random_state=rng).as_matrix()
    k = np.einsum('nqi,nij->nqj', np.pad(q, ((0,0),(0,0),(0,1))), r)
    # CTF-like radial phase transfer with varying defocus; not a microscope
    # simulator. Real microscope CTFs are reserved for the next experiment.
    radius2 = np.sum(q*q, axis=-1)
    transfer = -np.sin(rng.uniform(0.15, 0.45, (n,1))*radius2+0.08)
    return q, k, transfer


def posterior_methods(a, ell, B):
    precision = a.T @ a + (len(ell)/B**2)*np.eye(len(ell))
    v = cho_solve(cho_factor(precision), ell)
    w = a @ v
    z = norm.ppf(.975)
    return {
        'full_gaussian_posterior': (w, z*np.sqrt(ell@v)),
        'diagonal_variational_posterior': (w, z*np.sqrt(np.sum(ell**2/np.diag(precision)))),
        'ridge_sampling_covariance': (w, z*np.linalg.norm(w)),
    }


def make_target(centers, sigma, preferred):
    # Two fixed local averages represented by Gaussian point evaluations in
    # this finite model. Direction intentionally probes the missing-view axis.
    points = np.array([[.025, .015, .06], [.025, .015, -.06]])
    if not preferred:
        points = np.array([[.045, .025, .025], [-.045, -.025, -.025]])
    ell = np.diff(density_functionals(points, centers, sigma), axis=0)[0]
    return ell/np.linalg.norm(ell)


def run_case(seed, n, preferred, angle_deg, maxiter):
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    centers = rng.normal(size=(20,3))*2.2
    centers[centers[:,2]<0] *= -1
    sigma = 1.0
    p = 2*len(centers)
    pilot = rng.normal(size=p)
    pilot /= np.linalg.norm(pilot)
    B = .25
    ell = make_target(centers, sigma, preferred)
    q,k,transfer = geometry(rng,n,preferred)
    noise = .20
    angle = np.deg2rad(angle_deg)
    shift = 0 if angle_deg == 0 else .05
    box = 24
    ab = image_design(k,transfer,centers,sigma,noise)
    a = ab.reshape(-1,p)
    j = pose_jacobian(k,q,transfer,centers,sigma,pilot,box,angle,shift,noise)
    gamma,L,H = pose_remainder_bounds(k,q,transfer,centers,sigma,pilot,
                                     box,angle,shift,B,noise)
    derivatives=pose_design_derivatives(k,q,transfer,centers,sigma,box,angle,shift,noise)
    cross_design=derivatives.reshape(n,2*len(q[0]),-1)
    gamma2,_,_=pose_remainder_bounds(k,q,transfer,centers,sigma,pilot,
                                    box,angle,shift,B,noise,structured=True)
    methods = posterior_methods(a,ell,B)
    certificates = {}
    for name,ji,eta,gi,extra in [
        ('bias_aware_fixed_pose',np.zeros_like(j),0,0,[]),
        ('linearized_pose_only',j,1,0,[]),
        ('bounded_pose_remainder',j,1,gamma,[]),
        ('structured_pose_density',j,1,gamma2,[(cross_design,B)]),
    ]:
        begin = time.perf_counter()
        cert = optimize_certificate(a,ji,ell,B,eta,gi,maxiter=maxiter,rtol=2e-3,extra_nuisance_groups=extra)
        methods[name] = (cert.weights,cert.half_width)
        certificates[name] = {
            'noise_sd':cert.noise_sd,'density_bias':cert.reconstruction_bias,
            'pose_bias':cert.nuisance_bias,'remainder_bias':cert.remainder_bias,
            'objective':cert.objective,'dual_lower_bound':cert.dual_lower_bound,
            'relative_gap':(cert.objective-cert.dual_lower_bound)/max(cert.objective,1e-15),
            'converged':cert.converged,'iterations':cert.iterations,
            'nuisance_group_biases':cert.nuisance_group_biases,
            'seconds':time.perf_counter()-begin,
        }
    # A no-data interval is a mandatory width comparator, and has coverage one
    # for this known bounded class. Any useful robust interval should improve it.
    methods['density_ball_without_data'] = (np.zeros(a.shape[0]),B*np.linalg.norm(ell))
    records=[]
    baseline_w=methods['full_gaussian_posterior'][0]
    delta_direction=a.T@baseline_w-ell
    delta_direction/=np.linalg.norm(delta_direction)
    for stress in ['random_truth_pose','ridge_bias_boundary','joint_adversary','bounds_violated']:
        if stress=='random_truth_pose':
            delta=rng.normal(size=p);delta*=B/np.linalg.norm(delta)
            u=rng.normal(size=(n,5));u/=np.linalg.norm(u,axis=1,keepdims=True)
        else:
            delta=B*delta_direction
            u=np.einsum('nmq,nm->nq',j,baseline_w.reshape(n,-1))
            u/=np.maximum(np.linalg.norm(u,axis=1,keepdims=True),1e-300)
            if stress=='ridge_bias_boundary':
                u*=0
            if stress=='bounds_violated':
                delta*=2;u*=3
        expected=realify(perturbed_images(k,q,transfer,centers,sigma,pilot+delta,
                                         box,u,angle,shift,noise)).ravel()-a@pilot
        actual_remainder=(expected-a@delta).reshape(n,-1)-np.einsum('nmq,nq->nm',j,u)
        remainder_ratio=float(np.max(np.divide(np.linalg.norm(actual_remainder,axis=1),gamma,
                                  out=np.zeros(n),where=gamma>0)))
        truth=float(ell@delta)
        for name,(w,half) in methods.items():
            bias=float(w@expected-truth)
            sd=float(np.linalg.norm(w))
            coverage=float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd)) if sd else float(abs(bias)<=half+1e-12)
            records.append({'stress':stress,'method':name,'half_width':float(half),
                            'bias':bias,'noise_sd':sd,'exact_noise_coverage':coverage,
                            'width_fraction_of_no_data':float(half/B),
                            'max_actual_remainder_bound_ratio':remainder_ratio})
    return {'seed':seed,'particles':n,'frequencies':len(q[0]),'coefficients':p,
            'geometry':'preferred' if preferred else 'uniform','pose_degrees':angle_deg,
            'shift_pixels':shift,'density_radius':B,'known_noise_std':noise,
            'median_remainder_bound':float(np.median(gamma)),
            'median_density_pose_cross_bound':float(np.median(B*L)),
            'median_pilot_curvature_bound':float(np.median(H/2)),
            'seconds':time.perf_counter()-started,'certificates':certificates,'results':records}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--particles',type=int,default=96)
    parser.add_argument('--maxiter',type=int,default=120)
    parser.add_argument('--seed',type=int,default=20260929)
    parser.add_argument('--output',default='results/uncertainty/development/physical-pilot.json')
    args=parser.parse_args()
    destination=ROOT/args.output;destination.parent.mkdir(parents=True,exist_ok=True)
    output={'stage':'development: exact dictionary, known simulation bounds; not practical calibration',
            'python':platform.python_version(),'platform':platform.platform(),'cases':[]}
    for preferred in [False,True]:
        for angle in [0.,.5,2.]:
            result=run_case(args.seed,args.particles,preferred,angle,args.maxiter)
            output['cases'].append(result)
            destination.write_text(json.dumps(output,indent=2)+'\n')
            print(result['geometry'],angle,'seconds',round(result['seconds'],2),flush=True)
            for row in result['results']:
                if row['stress']=='joint_adversary':
                    print(row['method'],'coverage',round(row['exact_noise_coverage'],4),
                          'relative_width',round(row['width_fraction_of_no_data'],3),flush=True)


if __name__=='__main__':
    main()
