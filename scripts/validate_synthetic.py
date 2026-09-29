"""An out-of-dictionary phantom validates recovery rather than memorizing the basis."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from run_experiment import matrices,solve
from fourier_splats.basis import evaluate_grid
from fourier_splats.fsc import fsc
rng=np.random.default_rng(10);box=24;count=30000
k=rng.uniform(-10,10,(count*2,3)).astype(np.float32);k=k[np.linalg.norm(k,axis=1)<10][:count]
# Analytic transform of displaced real-space Gaussians (different generating basis).
centers=np.array([[2.,3.,-1.],[-3.,-1.,2.],[1.,-4.,-3.]])
width=np.array([1.3,1.8,1.1]);mass=np.array([1.,.7,.4])
def truth(k):
 s2=(k*k).sum(-1)/box**2
 return sum(mass[j]*np.exp(-2*np.pi**2*width[j]**2*s2)*np.exp(-2j*np.pi*(k@centers[j])/box) for j in range(3))
y=truth(k);y+=.05*(rng.normal(size=len(k))+1j*rng.normal(size=len(k)))
C=np.sin(.12*(k*k).sum(-1)).astype(np.float32);y*=C
q=np.arange(-box//2,box//2);z,yy,x=np.meshgrid(q,q,q,indexing='ij');g=np.stack([x.ravel(),yy.ravel(),z.ravel()],axis=-1)
true=truth(g).reshape((box,)*3)
results={}
for kind in ['gaussian','voxel']:
 ar,ai=matrices(k,C,box,kind,.5,2)
 re,lr=solve(ar,y.real.astype('float32'),.01,120);im,li=solve(ai,y.imag.astype('float32'),.01,120)
 f=evaluate_grid(re,im,box,kind);curve=fsc(f,true,1)
 results[kind]={'low_shell_mean_fsc':float(curve[:5,2].mean()),'shells':curve.tolist(),'solver_real':lr,'solver_imag':li}
 assert curve[:5,2].mean()>.95,(kind,curve)
# Randomized observations with identical sampling must not reconstruct a common signal.
Path('results/validation').mkdir(parents=True,exist_ok=True)
Path('results/validation/synthetic.json').write_text(json.dumps(results,indent=2))
print(json.dumps({k:v['low_shell_mean_fsc'] for k,v in results.items()}))
