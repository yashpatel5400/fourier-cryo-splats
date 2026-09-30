import numpy as np
from fourier_splats.uq_dual_mixture import closest_support_mixture


def test_support_mixture_agrees_with_exact_segment_projection():
    rng = np.random.default_rng(609901)
    for _ in range(8):
        directions = rng.normal(size=(2, 11)); target = rng.normal(size=11)
        delta = directions[1]-directions[0]
        t = np.clip((target-directions[0])@delta/(delta@delta), 0., 1.)
        fit = closest_support_mixture(directions, target, initial=[.3, .7])
        np.testing.assert_allclose(fit['mixture'], [1-t, t], atol=2e-8)
        np.testing.assert_allclose(fit['distance'], np.linalg.norm(target-directions[0]-t*delta), rtol=1e-12)


def test_redundant_modes_preserve_feasibility_and_never_worsen_support():
    directions = np.array([[1., 0., 0.], [0., 1., 0.], [1., 0., 0.], [0., 0., 1.]])
    fit = closest_support_mixture(directions, np.full(3, 1/3), initial=[1., 0., 0., 0.])
    assert fit['distance'] < 1e-7
    assert fit['distance'] <= fit['initial_distance']
    assert np.all(fit['mixture'] >= 0)
    np.testing.assert_allclose(fit['mixture'].sum(), 1.)


def test_joint_spectral_hull_and_pair_balls_matches_independent_conic_distance():
    import cvxpy as cp
    from fourier_splats.uq_dual_mixture import closest_support_and_pair_balls
    rng=np.random.default_rng(609971);n,nq,k=2,3,5;m=2*n*nq
    directions=rng.normal(size=(k,m));directions[-1]=0.;target=rng.normal(size=m);radii=.2+rng.random((n,nq))
    fit=closest_support_and_pair_balls(directions,target,radii)
    beta=cp.Variable(k);g=cp.Variable((n,2*nq))
    constraints=[beta>=0,cp.sum(beta)==1]
    for i in range(n):constraints.append(cp.norm(cp.vstack([g[i,:nq],g[i,nq:]]),axis=0)<=radii[i])
    problem=cp.Problem(cp.Minimize(cp.norm(target-beta@directions-cp.reshape(g,(m,),order='C'))),constraints)
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10)
    np.testing.assert_allclose(fit['distance'],problem.value,rtol=1e-6,atol=1e-7)
    np.testing.assert_allclose(fit['mixture'].sum(),1.,atol=1e-14)
    cubic=fit['cubic_support'].reshape(n,2*nq)
    assert np.all(np.hypot(cubic[:,:nq],cubic[:,nq:])<=radii+1e-14)
    np.testing.assert_allclose(fit['distance'],np.linalg.norm(target-fit['pose_support']-fit['cubic_support']),rtol=1e-12)
