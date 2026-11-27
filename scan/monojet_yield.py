"""Monojet-style yield of the neutral pair: charged quirk decayed and the pair leaves the detector (quirksim flag) AND generator-level
MET > cut, leading jet pT > 100 GeV, |eta_j1| < 2.4 (generic selection, no Delta-phi/jet-multiplicity/lepton-veto, 100% trigger above the cut).
usage: monojet_yield.py <ms dir> <out.csv>"""
import sys, numpy as np, pandas as pd
sys.argv += [] ; base, out = sys.argv[1], sys.argv[2]
MASSES = [300, 500, 750, 1000, 1500, 2000]; LAMS = [30, 100, 300, 1000]
KW2 = {'doublet': 1., 'triplet': 2.}; SPECIES = {'doublet': 1., 'triplet': 2.}
import re
def sigma_fb(f): return float(re.search(r'sigma = ([\d.eE+-]+) pb', open(f).read()).group(1)) * 1e3
rows = []
for M in MASSES:
    rw = np.load(f'{base}/M{M}/dy_rw.npz')
    for rep in ('doublet', 'triplet'):
        for L in LAMS:
            r = dict(M_GeV=M, rep=rep, Lambda_eV=L)
            for cut in (250., 400.):
                tot = 0.
                for p, fn in (('dy', 'dyf'), ('w', 'wf')):
                    ev = np.load(f'{base}/csv/{rep}_{p}_M{M}_ev.npz')[f'L{L}']; ids = ev[:, 0].astype(int); cover = ids.max() + 1
                    jm = np.load(f'{base}/M{M}/jm_{fn}.npz'); idx = {e: i for i, e in enumerate(jm['evnum'])}
                    k = np.array([idx[i] for i in ids]); met = jm['met'][k]; j1 = jm['jets'][k, 0]
                    sel = (met > cut) & (j1[:, 0] > 100.) & (np.abs(j1[:, 1]) < 2.4)
                    if p == 'dy':
                        wmap = dict(zip(rw['evnum'], rw[rep])); w = SPECIES[rep] * np.array([wmap[i] for i in ids])
                    else: w = np.full(len(ids), KW2[rep] * SPECIES[rep])
                    tot += sigma_fb(f'{base}/M{M}/{p}.log') / cover * (w * ev[:, 2] * sel).sum()
                r[f'sigma_monojet_MET{int(cut)}_fb'] = tot
            rows.append(r)
d = pd.DataFrame(rows); d.to_csv(out, index=False)
print(d[d.Lambda_eV == 100].round(5).to_string(index=False))
