"""Full-grid pose derivative audits through Fourier-moment adjoints.

Only twenty frequency-weighted Fourier moment fields per particle are needed
for first and second pose derivatives, independent of fitting dictionary size.
They can be evaluated by direct summation or batched type-1 NUFFTs.
The returned norms audit fixed weights; this module does not optimize weights.
"""
import numpy as np
import finufft


def pose_derivative_adjoints(operator,q,weights,rotation_radius,shift_radius,particle,backend='auto'):
    """Return D_a' w and E_ab' w in the supported voxel coordinates."""
    op=operator.op;i=particle;xyz=operator.xyz;nq=op.q
    wi=np.asarray(weights).reshape(op.n,2*nq)[i]
    frequencies=np.concatenate([op.k[i],np.asarray(q)[i]],axis=1)
    pairs=[(a,b) for a in range(5) for b in range(a,5)]
    moments=np.concatenate([frequencies.T,np.stack([frequencies[:,a]*frequencies[:,b] for a,b in pairs])])
    coef=(wi[:nq]+1j*wi[nq:])*op.transfer[i]
    coefficients=np.ascontiguousarray(moments*coef[None],dtype=complex)
    if backend=='auto':backend='direct' if nq<=128 else 'nufft'
    if backend=='direct':
        transform=np.empty((len(xyz),20),complex)
        for start in range(0,len(xyz),8192):
            phase=np.exp(2j*np.pi*(op.k[i]@xyz[start:start+8192].T)/op.box)
            transform[start:start+len(phase.T)]=(coefficients@phase).T
    elif backend=='nufft':
        coordinates=2*np.pi*op.k[i]/op.box
        transform=finufft.nufft3d1(*[np.ascontiguousarray(coordinates[:,a]) for a in [2,1,0]],
            coefficients,(op.box,)*3,isign=1,eps=op.eps,nthreads=op.nthreads)
        transform=transform.reshape(20,-1)[:,operator.active].T
    else:raise ValueError('Backend must be auto, direct or nufft')
    first=transform[:,:5]
    second=np.empty((len(xyz),5,5),complex)
    for idx,(a,b) in enumerate(pairs):second[:,a,b]=second[:,b,a]=transform[:,5+idx]
    scale=2*np.pi/op.box
    # Maps the five frequency components to the five pose-phase derivatives.
    linear=np.zeros((len(xyz),5,5))
    x,y,z=xyz.T
    linear[:,0,1]=-scale*rotation_radius*z;linear[:,0,2]=scale*rotation_radius*y
    linear[:,1,0]=scale*rotation_radius*z;linear[:,1,2]=-scale*rotation_radius*x
    linear[:,2,0]=-scale*rotation_radius*y;linear[:,2,1]=scale*rotation_radius*x
    linear[:,3,3]=scale*shift_radius;linear[:,4,4]=scale*shift_radius
    d1=(1j*np.einsum('par,pr->pa',linear,first)).real
    d2=-np.einsum('par,prs,pbs->pab',linear,second,linear).real
    dot=np.sum(xyz*first[:,:3],axis=1)
    for a in range(3):
        for b in range(3):
            curvature=scale*rotation_radius**2*(.5*(xyz[:,b]*first[:,a]+xyz[:,a]*first[:,b])-(a==b)*dot)
            d2[:,a,b]+=(1j*curvature).real
    return d1,d2


def audit_pose_weights(operator,pilot,q,weights,rotation_radius,shift_radius,density_radius,order=2,backend='auto'):
    """Evaluate all structured nuisance norms and uniform remainder on this grid."""
    if order not in [1,2]:raise ValueError('Taylor order must be one or two')
    if not 0<=rotation_radius<=np.pi or shift_radius<0 or density_radius<0:raise ValueError('Invalid radii')
    op=operator.op;pilot=np.asarray(pilot);q=np.asarray(q);w=np.asarray(weights).reshape(op.n,2*op.q)
    a=rotation_radius;s=shift_radius;B=density_radius;c=2*np.pi/op.box
    nx=np.linalg.norm(operator.xyz,axis=1);rows=[]
    for i in range(op.n):
        d1,d2=pose_derivative_adjoints(operator,q,w,a,s,i,backend)
        linear_pilot=float(np.linalg.norm(pilot@d1));linear_density=float(B*np.linalg.norm(d1))
        quadratic_pilot=float(.5*np.linalg.norm(np.einsum('p,pab->ab',pilot,d2))) if order==2 else 0.
        quadratic_density=float(.5*B*np.linalg.norm(d2)) if order==2 else 0.
        kr=np.linalg.norm(op.k[i],axis=1)
        speed=c*(a*kr[:,None]*nx[None,:]+s*np.linalg.norm(q[i],axis=1)[:,None])
        acceleration=c*a*a*kr[:,None]*nx[None,:]
        if order==1:
            column=np.abs(op.transfer[i,:,None])*(speed**2+acceleration);factor=.5
        else:
            jerk=c*a**3*kr[:,None]*nx[None,:]
            column=np.abs(op.transfer[i,:,None])*(speed**3+3*speed*acceleration+jerk);factor=1/6
        gamma=factor*(np.linalg.norm(column@np.abs(pilot))+B*np.linalg.norm(column))
        rows.append([linear_pilot,linear_density,quadratic_pilot,quadratic_density,float(gamma*np.linalg.norm(w[i]))])
    names=['linear_pilot','linear_density','quadratic_pilot','quadratic_density','remainder']
    return {'per_particle':np.asarray(rows),'totals':dict(zip(names,map(float,np.sum(rows,axis=0))))}
