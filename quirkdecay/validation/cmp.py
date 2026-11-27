import sys, uproot, numpy as np
cols=['genEvt','maxSep','q_xe','q_ye','q_ze','q_te','q_ee']
def load(f):
    a=uproot.open(f)['events'].arrays(cols,library='np'); i=np.argsort(a['genEvt']); return {k:np.array([a[k][j] for j in i],float) for k in cols}
A,B=load(sys.argv[1]),load(sys.argv[2])
for k in cols[1:]:
    d=np.abs(A[k]-B[k]); rel=d/(np.abs(B[k])+1e-9)
    print(f'{sys.argv[1]} vs {sys.argv[2]} {k}: max|diff| {d.max():.3g} median rel {np.median(rel):.2e} 95% {np.percentile(rel,95):.2e} max {rel.max():.2e}')
