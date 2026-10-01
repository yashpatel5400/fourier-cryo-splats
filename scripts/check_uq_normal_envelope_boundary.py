#!/usr/bin/env python3
"""Independent scalar checks of two sharp Gaussian-envelope boundaries."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/uncertainty/development/audit-regressions/sharp-normal-envelope.json'
    if out.exists():raise RuntimeError('Preserve existing numerical check')
    cases=[]
    for z in [1.,1.5,1.7,np.sqrt(3),norm.isf(.025),norm.isf(.05/24)]:
        t=.025;q=z*np.sqrt(1+t*t)
        observed=norm.cdf(q-t)+norm.cdf(q+t)-2*norm.cdf(z)
        coefficient=norm.pdf(z)*z*(z*z-3)/6
        if z<np.sqrt(3) and not observed<0:raise AssertionError('Expected small-bias counterexample')
        if z>np.sqrt(3) and not observed>0:raise AssertionError('Expected positive envelope difference')
        cases.append(dict(z=float(z),t=t,coverage_difference=float(observed),
            t4_coefficient=float(coefficient),difference_over_t4=float(observed/t**4)))
    z=norm.isf(.025);b=1/np.sqrt(z*z-1)
    grid=np.linspace(0,3,10001)
    empirical_min=float(np.min(2*norm.cdf(z*np.sqrt(1+grid**2)-grid)-1))
    exact=float(2*norm.cdf(np.sqrt(z*z-1))-1)
    if abs(empirical_min-exact)>1e-8:raise AssertionError('Grid disagrees with exact adaptive-bias minimum')
    d=dict(complete=True,scope='Scalar numerical support for an elementary derivation, not scientific novelty or data-adaptive pose validity.',
        fixed_bias_cases=cases,nominal=.95,adaptive_bias=dict(z=float(z),worst_b=float(b),
            exact_minimum_coverage=exact,grid_minimum_coverage=empirical_min,
            valid_rms_constant=float(np.sqrt(1+z*z)),relative_constant_increase=float(np.sqrt(1+z*z)/z-1)))
    out.write_text(json.dumps(d,indent=2)+'\n');print(d['adaptive_bias'])


if __name__=='__main__':main()
