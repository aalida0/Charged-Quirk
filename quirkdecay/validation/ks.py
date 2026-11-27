import uproot, numpy as np
from scipy import stats
cols=['genEvt','maxSep','q_xe','q_ye','q_ze','q_te','q_ee']
def load(f):
    a=uproot.open(f)['events'].arrays(cols,library='np'); i=np.argsort(a['genEvt']); return {k:np.array([a[k][j] for j in i],float) for k in cols}
D={t:load(t+'.root') for t in ['vA','vB','vC','vI2','vJ100i','vJ1000i']}
print('KS p-values and mean shift (in units of stat. error) for A(decay,Δ=2MeV) vs B(neutral from start); and identity-vs-no-decay control')
for a,b in [('vA','vB'),('vI2','vC'),('vJ100i','vC'),('vJ1000i','vC')]:
    for k in ['maxSep','q_ze','q_te','q_ee']:
        x,y=D[a][k].ravel(),D[b][k].ravel(); ks=stats.ks_2samp(x,y)
        d=x-y; z=d.mean()/(d.std(ddof=1)/np.sqrt(len(d)))
        print(f'{a:7s} vs {b:4s} {k:7s} KS p={ks.pvalue:.2f}  paired mean diff {d.mean():+.3g} ({z:+.1f} sigma)  mean {y.mean():.4g}')
