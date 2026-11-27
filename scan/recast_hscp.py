"""Survival-weighted HSCP cross section of decaying charged quirks.
usage: recast_hscp.py <ms dir> <out.csv>
For every charged quirk of the mass-scan samples (generator kinematics from the quirksim ntuple, DY re-weighted + W*), the probability to
reach radius R before decaying is exp(-L_R/(beta*gamma*c*tau)), with L_R the straight-line path from the production point to the
cylinder r = R (|z| < Zmax); string slowing and magnetic bending are neglected (check: scan/ loss_check shows the median kinetic-energy loss
to the string before the decay is a few percent).  sigma_HSCP = sigma * P(at least one charged quirk reaches R).  Output: P and sigma_HSCP
per mass, representation, R and, for the proper c*tau actually predicted and for a c*tau scan.
To obtain a recast limit, compare sigma_HSCP(M) with the published upper limit sigma_UL(M) of the track-based HSCP search (not included)."""
import sys, glob, re, csv, numpy as np, uproot
base, out = sys.argv[1], sys.argv[2]
MASSES = [300, 500, 750, 1000, 1500, 2000]; KW2 = {'doublet': 1., 'triplet': 2.}
RADII = {'tracker': (1160., 2800.), 'ecal': (1290., 3150.), 'muon': (4200., 6900.)}
CTAUS = [1., 5., 20., 56., 200., 1000., 5000., 20000., 1e5]   # mm
def sigma_fb(f): return float(re.search(r'sigma = ([\d.eE+-]+) pb', open(f).read()).group(1)) * 1e3
def Lexit(x0, y0, z0, px, py, pz, R, Z):
    pt2 = px * px + py * py; pp = np.sqrt(pt2 + pz * pz); ux, uy, uz = px / pp, py / pp, pz / pp
    a = ux * ux + uy * uy; b = 2 * (x0 * ux + y0 * uy); c = x0 * x0 + y0 * y0 - R * R
    s_r = (-b + np.sqrt(max(b * b - 4 * a * c, 0.))) / (2 * a) if a > 0 else np.inf
    s_z = ((Z - z0) / uz if uz > 0 else (-Z - z0) / uz) if uz != 0 else np.inf
    return min(s_r, s_z)
rows = []
for M in MASSES:
    rw = np.load(f'{base}/M{M}/dy_rw.npz'); sg = {p: sigma_fb(f'{base}/M{M}/{p}.log') for p in ('dy', 'w')}
    for rep in ('doublet', 'triplet'):
        ct_phys = float(re.search(rf'M={M} {rep} .*ctau=([\d.]+) mm', open(f'{base}/ms_params.txt').read()).group(1))
        data = {}
        for p in ('dy', 'w'):
            a = uproot.open(f'{base}/out/{rep}_{p}_M{M}/L100_c{ct_phys:.4f}.root')['events'].arrays(
                ['genEvt', 'q_pdg', 'q_px0', 'q_py0', 'q_pz0', 'q_x0', 'q_y0', 'q_z0'], library='np')
            ids = a['genEvt'].astype(int); cover = ids.max() + 1
            wmap = dict(zip(rw['evnum'], rw[rep])) if p == 'dy' else None
            w = np.array([wmap[i] for i in ids]) if p == 'dy' else np.full(len(ids), KW2[rep])
            bg, path = [], {k: [] for k in RADII}
            for i in range(len(ids)):
                bgs, ls = [], {k: [] for k in RADII}
                for s in (0, 1):
                    if abs(int(a['q_pdg'][i][s])) != 17: continue
                    px, py, pz = a['q_px0'][i][s], a['q_py0'][i][s], a['q_pz0'][i][s]
                    bgs.append(np.sqrt(px * px + py * py + pz * pz) / (M * 1e3))
                    for k, (R, Z) in RADII.items(): ls[k].append(Lexit(a['q_x0'][i][s], a['q_y0'][i][s], a['q_z0'][i][s], px, py, pz, R, Z))
                bg.append(bgs)
                for k in RADII: path[k].append(ls[k])
            data[p] = (w, bg, path, sg[p] / cover)
        for ct in [ct_phys] + CTAUS:
            for k in RADII:
                sig = 0.; stot = 0.
                for p in ('dy', 'w'):
                    w, bg, path, sc = data[p]
                    P = np.array([1 - np.prod([1 - np.exp(-l / (b * ct)) for l, b in zip(path[k][i], bg[i])]) for i in range(len(w))])
                    sig += sc * (w * P).sum(); stot += sc * w.sum()
                rows.append(dict(M_GeV=M, rep=rep, ctau_mm=ct, physical=int(ct == ct_phys), radius=k, sigma_fb=stot, sigma_HSCP_fb=sig, P_mean=sig / stot))
    print('done', M, flush=True)
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows:
    if r['physical']: print(f"{r['rep']:7s} M={r['M_GeV']:5d} ct={r['ctau_mm']:.1f}mm {r['radius']:8s} P={r['P_mean']:.3g} sigma_HSCP={r['sigma_HSCP_fb']:.3g} fb")
