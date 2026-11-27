#!/usr/bin/env python3
"""Transport closure of the UFO-based recoil sample against the Pythia-template samples (same M, Lambda, decay table), per class (DY / W*).
Per charged-quirk endpoint that decayed in the tracker (|q_pdg|=17): fraction with >= 3 / >= 4 distinct tracker layers (Geant4 sensitive-volume hits).
usage: closure_transport.py <ms dir> <M> <Lambda list comma>"""
import sys, glob, numpy as np, uproot
base, M, Ls = sys.argv[1], int(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')]
TRK = [2, 4, 5, 6, 7, 8]
cls = {}; cur = None
for l in open(f'{base}/M{M}/rec_ufo.hepmc3'):
    if l.startswith('E '): cur = int(l.split()[1]); cls[cur] = 'dy'
    elif l.startswith('P ') and l.split()[-1] == '1' and abs(int(l.split()[3])) == 18: cls[cur] = 'w'
def endpoints(fn, select=None):
    a = uproot.open(fn)['events'].arrays(['genEvt', 'q_pdg', 'dec_side', 'h_q', 'h_det', 'h_layer'], library='np')
    n_ev = n_dec = n3 = n4 = 0
    for i, g in enumerate(a['genEvt'].astype(int)):
        if select is not None and cls.get(g) != select: continue
        n_ev += 1
        for s in {int(v) for v in a['dec_side'][i]}:
            if abs(int(a['q_pdg'][i][s])) != 17: continue
            n_dec += 1
            sel = (a['h_q'][i] == s) & np.isin(a['h_det'][i], TRK)
            nl = len({(int(d), int(l)) for d, l in zip(a['h_det'][i][sel], a['h_layer'][i][sel])})
            n3 += nl >= 3; n4 += nl >= 4
    return n_ev, n_dec, n3, n4
for L in Ls:
    print(f'--- M={M} Lambda={L} eV')
    for c in ('dy', 'w'):
        ref = endpoints(glob.glob(f'{base}/out3/triplet_{c}_M{M}/L{L}_c*.root')[0])
        uf = endpoints(glob.glob(f'{base}/out3/triplet_ufo_M{M}/L{L}_c*.root')[0], select=c)
        for nm, r in (('Pythia template', ref), ('UFO + Pythia', uf)):
            ne, nd, n3, n4 = r
            print(f'  {c:2s} {nm:15s}: events {ne:5d}, charged endpoints decayed {nd:5d} ({nd/max(ne,1):.3f}/event), >=3 layers {n3/max(nd,1):.3f}+-{np.sqrt(n3/max(nd,1)*(1-n3/max(nd,1))/max(nd,1)):.3f}, >=4 layers {n4/max(nd,1):.3f}')
