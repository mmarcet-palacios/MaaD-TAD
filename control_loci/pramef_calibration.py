"""pramef_calibration.py -- homologous-pair contact enrichment in PRAMEF at the calibrated homology weight (Fig 9C, right).
Run in MaaD-TAD\\control_loci after sweep_calibrate.py:   python pramef_calibration.py
Folds PRAMEF_seqsim_like.csv with csynth_like.fold at k = 0 and k = 0.007 (homologous pairs only, five starts) and
measures the same enrichment (median distance^-3 ratio, homologous vs same-separation control) in the CSynth path,
the refolded path and the final (tethered) path from ..\\data. Writes pramef_cal.json."""
import numpy as np, json, sys, os
HERE=os.path.dirname(os.path.abspath(__file__)); DATA=os.path.join(HERE,'..','data')
from csynth_like import fold
from multiprocessing import Pool
n=325
S=np.loadtxt(os.path.join(HERE,'PRAMEF_seqsim_like.csv'),delimiter=',',skiprows=1,converters={0:lambda s:0},usecols=range(1,n+1))[:n,:n]
iu=np.triu_indices(n,3); sep=iu[1]-iu[0]; sv=S[iu]; hom=sv>1.1; rng=np.random.default_rng(0); hk=np.where(hom)[0]
ctrl=np.array([rng.choice(np.where((sv<=0.6)&(sep==sep[i]))[0]) for i in hk])
def enr(X):
    D=np.linalg.norm(X[:,None]-X[None],axis=-1); P=D[iu]**-3.0; return float(np.median(P[hk]/P[ctrl]))
res={'n_hom_pairs':int(len(hk))}
for name in ['gap-filled','guide_refolded','guide_tethered']:
    X=np.loadtxt(os.path.join(DATA,f'{name}.xyz'))[:n,1:]; res[name]=enr(X)
def job(a):
    k,seed=a; M=S if k>0 else np.full((n,n),0.4)
    X,_=fold(M,kc=(k if k>0 else 0.25),rep=0.3,seeds=1,rng0=100+seed,thresh=1.1); return k,seed,enr(X)
if __name__=='__main__':
  with Pool(2) as p: rows=p.map(job,[(k,s) for k in (0,0.007) for s in range(5)])
  for k in (0,0.007):
      e=[r[2] for r in rows if r[0]==k]; res[f'csynth_like_k{k}']=(float(np.mean(e)),float(np.std(e)))
  print(json.dumps(res,indent=1)); json.dump(res,open(os.path.join(HERE,'pramef_cal.json'),'w'),indent=1)
