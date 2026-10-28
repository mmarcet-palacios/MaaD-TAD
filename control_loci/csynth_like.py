"""CSynth-like polymer folding from a similarity/contact matrix.
Beads = 2 kb bins. Energy = backbone springs (rest length 1) + contact attraction w_ij*d_ij^2
+ soft pairwise repulsion rep/(d^2+0.25). Minimised with L-BFGS from several random-walk starts."""
import numpy as np
from scipy.optimize import minimize
def fold(C, kb=50.0, kc=1.0, rep=0.3, bg=None, seeds=4, iters=4000, rng0=0, thresh=None):
    n=len(C); C=np.array(C,float)
    if bg is None: bg=np.nanmedian(C[np.triu_indices(n,2)])
    W=np.nan_to_num(C,nan=bg)/100.0 - bg/100.0; W[W<0]=0; np.fill_diagonal(W,0)
    if thresh is not None: W[np.nan_to_num(C,nan=0)<thresh]=0   # attract only pairs called homologous
    for k in (1,-1): W[np.arange(max(0,-k),n-max(0,k)),np.arange(max(0,k),n+min(0,k))]=0   # neighbours: backbone only
    W=kc*W
    iu=np.triu_indices(n,1); Wu=W[iu]; far=(iu[1]-iu[0])>1
    def E(xf):
        X=xf.reshape(n,3); D=X[:,None,:]-X[None,:,:]; d2=(D**2).sum(-1)
        # backbone
        b=X[1:]-X[:-1]; bl=np.sqrt((b**2).sum(1))+1e-9
        e=kb*((bl-1)**2).sum(); g=np.zeros_like(X)
        gb=(2*kb*(bl-1)/bl)[:,None]*b; g[1:]+=gb; g[:-1]-=gb
        # contacts (harmonic to 0) + repulsion
        d2u=d2[iu]; e+= (Wu*d2u).sum()
        r=rep*far/(d2u+0.25); e+=r.sum()
        coef=2*Wu - 2*rep*far/(d2u+0.25)**2          # dE/d(d2) * 2
        G=np.zeros((n,n)); G[iu]=coef; G=G+G.T
        g+= (G[:,:,None]*D).sum(1)
        return e,g.ravel()
    best=None
    for s in range(seeds):
        rng=np.random.default_rng(rng0+s); X0=np.cumsum(rng.normal(size=(n,3))/np.sqrt(3),0)
        r=minimize(E,X0.ravel(),jac=True,method='L-BFGS-B',options=dict(maxiter=iters))
        if best is None or r.fun<best.fun: best=r
    return best.x.reshape(n,3), best.fun
