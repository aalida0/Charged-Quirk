"""Mass limit of the wino chain at tau=0.19 ns (the quirk-benchmark lifetime) from sigma_vis(M)=sum_ch sigma_NLO+NLL,ch(M) eff_ch(M) = 0.037 fb,
for the 'first' and 'intended' smearing conventions, against the official-contour value (649 GeV at ctau=56 mm). usage: wino019_limit.py <data/massscan>"""
import sys, numpy as np, pandas as pd
dd = sys.argv[1]
ev = pd.concat([pd.read_csv(f'{dd}/wino_validation019_{p}.csv') for p in ('chgchg', 'chgn1')])
xs = {p: pd.read_csv(f'{dd}/xsec_13TeV_pp2{p}_fb.csv').set_index('mass').xsec for p in ('chgchg', 'chgn1')}
grid = np.linspace(450, 850, 801)
for sm in ('first', 'intended'):
    vis = np.zeros_like(grid)
    for p in ('chgchg', 'chgn1'):
        e = ev[(ev.proc == p) & (ev.smear == sm)].sort_values('M'); Ms = e.M.values.astype(float)
        sg = np.exp(np.interp(grid, xs[p].index.values.astype(float), np.log(xs[p].values)))
        vis += sg * np.interp(grid, Ms, e.eff.values)
    k = np.where(vis > 0.037)[0]; print(f'{sm}: sigma_vis(600)={np.interp(600, grid, vis):.4f} fb; limit = {grid[k].max():.0f} GeV (official contour 649 GeV)')
