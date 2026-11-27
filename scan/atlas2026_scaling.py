"""Scaling estimate (NOT a recast) of the reach of the 2026 ATLAS disappearing-track search (arXiv:2603.08315, 137/fb) for helical (Lambda -> 0)
Dirac-multiplet quirks, from the benchmark limits quoted in the paper: wino-like (Delta m ~ 160 MeV, tau ~ 0.2 ns): 670 GeV observed (770 GeV expected);
higgsino-like (Delta m 280-350 MeV, tau 0.026-0.048 ns): 225 GeV observed (250 GeV expected).
Assumption: same efficiency as the benchmark, i.e. the visible rate scales with the production cross section only; the limit mass M* then solves
   mult * sigma_ref(M*) = sigma_ref(M_bench)        with mult = (rate of the quirk multiplet / rate of the benchmark per infracolor index) x N_ic.
Triplet: Dirac triplet = 2 x wino (per index). Doublet: Dirac doublet = higgsino-like; sigma_higgsino(M) = r x sigma_wino(M) with r = 0.41 at 300 GeV from the
paper's reference cross sections (higgsino 90.8+90.8+53 fb vs wino 387+190 fb) assumed mass independent (+-20%), which cancels in M* because both benchmark and quirk use the same shape.
sigma_wino(M): public NLO+NLL tables (data/massscan/xsec_13TeV_pp2chg*_fb.csv).   usage: atlas2026_scaling.py <data dir>"""
import sys, numpy as np, pandas as pd
dd = sys.argv[1]
xs = sum(pd.read_csv(f'{dd}/xsec_13TeV_pp2{p}_fb.csv').set_index('mass').xsec for p in ('chgchg', 'chgn1'))
M, S = xs.index.values.astype(float), xs.values
sig = lambda m: float(np.exp(np.interp(m, M, np.log(S))))
def mstar(mult, mb):
    f = lambda m: mult * sig(m) - sig(mb)
    grid = np.linspace(M.min(), M.max(), 4000); v = np.array([f(m) for m in grid]); k = np.where(v > 0)[0].max()
    return grid[k] if k < len(grid) - 1 else np.nan
print('wino benchmark 670 (obs) / 770 (exp) GeV, Dirac triplet = 2 winos per index:')
for nic in (1, 2, 3): print(f'  N_ic={nic}: M* = {mstar(2*nic, 670):.0f} (obs-scaled)  {mstar(2*nic, 770):.0f} (exp-scaled) GeV')
print('higgsino benchmark 225 (obs) / 250 (exp) GeV, Dirac doublet = 1 higgsino multiplet per index:')
for nic in (1, 2, 3): print(f'  N_ic={nic}: M* = {mstar(nic, 225):.0f} (obs-scaled)  {mstar(nic, 250):.0f} (exp-scaled) GeV')
