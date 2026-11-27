"""Filter-closure test for the pair-pT > 180 GeV pre-filter of the recoil-enriched samples (referee C7):
unfiltered (|eta| only) samples, generator-level jets/MET with the SAME smearing as atlas_dt_recast.py, trigger-weighted pre-selection.
Reports the trigger-weighted fraction of the passing events that has generator pair pT < PTCUT.
usage: filter_closure.py <ms dir> <ATLAS_data dir> <kfactor_check.csv> [PTCUT=180]"""
import sys, re, numpy as np, uproot
PT = float(sys.argv[4]) if len(sys.argv) > 4 else 180.
src = open(sys.argv[0].replace('filter_closure.py', 'atlas_dt_recast.py')).read()
head = src[:src.index('params = {}')]
sys.argv = [sys.argv[0]] + sys.argv[1:4] + ['/dev/null', 'intended']
exec(head)
for M in (300, 500, 750, 1000, 1500):
    for p in ('dy', 'w'):
        try:
            z = np.load(f'{base}/M{M}/{p}.npz'); c = [str(x) for x in z['cols']]; d = z['data']; jm = np.load(f'{base}/M{M}/jm2_{p}f.npz')
        except Exception as e: continue
        ev = {int(e): k for k, e in enumerate(d[:, 0])}; tot = lo = nlo = n = 0
        for k, e in enumerate(jm['evnum']):
            if int(e) not in ev: continue
            pr = d[ev[int(e)], c.index('pr_pt')]
            pt, eta, ph, met, mph = event_objects(jm, k); trig, kin = presel(pt, eta, ph, met, mph)
            if trig <= 0 or not kin: continue
            tot += trig; n += 1
            if pr < PT: lo += trig; nlo += 1
        print(f'M={M} {p}: events passing {n}, of which with pair pT<{PT:.0f}: {nlo}  -> trigger-weighted fraction lost by the filter = {lo / max(tot, 1e-9):.4f}')
