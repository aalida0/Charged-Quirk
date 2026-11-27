"""[v2, referee-v3 P0-B] Joint-acceptance pilot (referee item C4): trigger x kinematic pre-selection (ATLAS-SUSY-2018-19 map, generator-level jets/MET with the
parametrised smearing of atlas_dt_recast.py) x hit multiplicity of the decayed charged quirk (DISTINCT tracker layers, one hit per layer)
x measurable non-helicity, at the physical (loop-predicted) triplet lifetime, on the actual recoil-selected samples, NLO+NLL weights.
Per charged-quirk endpoint with >= 4 distinct layers a smeared circle fit (10 um pixel, 30 um strip, geometric chi2) is made. Three SEPARATE quantities are
recorded: chi2/ndf>3 and chi2/ndf>9 (non-circularity) and R<0.5 m (strong curvature, NOT a non-helicity criterion). A matched string-free control (Lambda=1 eV,
same samples, same code) gives the conditional false-positive rate of the chi2 flag at fixed >=4-layer acceptance. Weighted yields carry N_eff=(sum w)^2/sum w^2.
Three layers define a circle exactly and cannot be tested (>=3-layer column only). Lambda points are those with a finished simulation (.done).
Hits are Geant4 sensitive-volume visits of the CMS-like generic tracker (4+ pixel/strip detectors as in the card), NOT reconstructed tracks, no
pattern recognition, no fakes/backgrounds. Event-level: an event counts once if any decayed charged quirk qualifies.
usage: pilot_joint.py <ms dir> <ATLAS_data dir> <kfactor_check.csv> <out.csv>"""
import sys, os, re, csv, glob, numpy as np, uproot
src = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'atlas_dt_recast.py')).read()
head = src[:src.index('params = {}')]
sys.argv = [sys.argv[0]] + sys.argv[1:5] + ['intended']
exec(head)                                   # ATLAS maps, smearing, K(M), event_objects, presel, quirks, track_eff
OUT = sys.argv[4]
TRK = [2, 4, 5, 6, 7, 8]; PIX = {2, 4}
def cfit(x, y, s):
    w = 1. / s ** 2
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]; b = x ** 2 + y ** 2
    c, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], b * np.sqrt(w), rcond=None)
    xc, yc = c[0], c[1]; R = np.sqrt(max(c[2] + xc ** 2 + yc ** 2, 1e-12))
    for _ in range(5):
        dx, dy = x - xc, y - yc; rr = np.hypot(dx, dy) + 1e-12
        J = np.c_[-dx / rr, -dy / rr, -np.ones_like(x)] / s[:, None]
        step, *_ = np.linalg.lstsq(J, -(rr - R) / s, rcond=None); xc, yc, R = xc + step[0], yc + step[1], R + step[2]
    res = (np.hypot(x - xc, y - yc) - R) / s
    return float(np.sum(res ** 2)) / (len(x) - 3), abs(R)
rng2 = np.random.default_rng(77)
def layers(a, i, s):
    """first hit of each distinct tracker layer of quirk s, time ordered: arrays x, y, sigma"""
    sel = np.where((a['h_q'][i] == s) & np.isin(a['h_det'][i], TRK))[0]
    if len(sel) == 0: return None
    sel = sel[np.argsort(a['h_t0'][i][sel])]
    seen, keep = set(), []
    for k in sel:
        key = (int(a['h_det'][i][k]), int(a['h_layer'][i][k]))
        if key in seen: continue
        seen.add(key); keep.append(k)
    keep = np.array(keep)
    x = ((a['h_x0'][i] + a['h_x1'][i]) / 2)[keep]; y = ((a['h_y0'][i] + a['h_y1'][i]) / 2)[keep]
    sg = np.where(np.isin(a['h_det'][i][keep], list(PIX)), 0.010, 0.030)
    return x + rng2.normal(0, sg), y + rng2.normal(0, sg), sg
params = {}
for l in open(f'{base}/ms_params.txt'):
    m = re.match(r'M=(\d+) (\w+) dm=([\d.eE+-]+) GeV ctau=([\d.]+) mm', l); params[(int(m[1]), m[2])] = (float(m[3]), float(m[4]))
rows = []; rep = 'triplet'
MASSES = [int(m) for m in os.environ.get('PJ_M', '300,500,750,1000,1500').split(',')]
for M in MASSES:
    sg_, ng_ = {}, {}
    for p in ('dy', 'w'):
        t = open(f'{base}/M{M}/rec_{p}.log').read(); sg_[p] = float(re.search(r'sigma = ([\d.eE+-]+) pb', t).group(1)) * 1e3; ng_[p] = int(re.search(r'wrote (\d+) events', t).group(1))
    rw = np.load(f'{base}/M{M}/rec_dy_rw.npz'); jmz = {p: np.load(f'{base}/M{M}/jm3_{p}.npz') for p in ('dy', 'w')}
    Ls = sorted({int(re.search(r'/L(\d+)_c', f)[1]) for f in glob.glob(f'{base}/out3/{rep}_dy_M{M}/L*_c*.done')} &
                {int(re.search(r'/L(\d+)_c', f)[1]) for f in glob.glob(f'{base}/out3/{rep}_w_M{M}/L*_c*.done')})
    for L in Ls:
        keys = ('presel', 'L3', 'L4', 'C3', 'C9', 'R5', 'C3x2', 'C3_or_R5')
        S = {k: 0. for k in keys}; E2 = {k: 0. for k in keys}; N = {k: 0 for k in keys}; nep = 0; nep3 = 0; nep9 = 0; nepR = 0
        for p in ('dy', 'w'):
            fn = glob.glob(f'{base}/out3/{rep}_{p}_M{M}/L{L}_c*.root')[0]
            a = uproot.open(fn)['events'].arrays(['genEvt', 'q_pdg', 'dec_side', 'h_q', 'h_det', 'h_layer', 'h_x0', 'h_y0', 'h_x1', 'h_y1', 'h_t0'], library='np')
            ids = a['genEvt'].astype(int); scale = sg_[p] / ng_[p] * K(M, p)
            if p == 'dy': wmap = dict(zip(rw['evnum'], rw[rep])); w = SPECIES[rep] * np.array([wmap[i] for i in ids])
            else: w = np.full(len(ids), KW2[rep] * SPECIES[rep])
            jm = jmz[p]; idx = {int(e): k for k, e in enumerate(jm['evnum'])}
            for i, g in enumerate(ids):
                pt, eta, ph, met, mph = event_objects(jm, idx[g]); trig, kin = presel(pt, eta, ph, met, mph)
                if trig <= 0. or not kin: continue
                wt = scale * w[i] * trig
                flags = {k: False for k in keys}; flags['presel'] = True; n_c3 = 0
                for s in {int(v) for v in a['dec_side'][i]}:
                    if abs(int(a['q_pdg'][i][s])) != 17: continue
                    ly = layers(a, i, s)
                    if ly is None: continue
                    n = len(ly[0])
                    if n >= 3: flags['L3'] = True
                    if n >= 4:
                        flags['L4'] = True; nep += 1
                        c2, R = cfit(*ly)
                        f3, f9, fR = c2 > 3., c2 > 9., R < 500.
                        nep3 += f3; nep9 += f9; nepR += fR
                        if f3: flags['C3'] = True; n_c3 += 1
                        if f9: flags['C9'] = True
                        if fR: flags['R5'] = True
                        if f3 or fR: flags['C3_or_R5'] = True
                if n_c3 >= 2: flags['C3x2'] = True
                for k in keys:
                    if flags[k]: S[k] += wt; E2[k] += wt * wt; N[k] += 1
        dm, ct = params[(M, rep)]
        row = dict(M_GeV=M, Lambda_eV=L, ctau_mm=ct, nEndpoints_L4=nep, nEndpoints_chi2gt3=nep3, nEndpoints_chi2gt9=nep9, nEndpoints_Rlt500=nepR)
        for k in keys:
            row[f'sigma_{k}_fb'] = S[k]; row[f'err_{k}_fb'] = E2[k] ** .5; row[f'Nraw_{k}'] = N[k]
            row[f'Neff_{k}'] = S[k] ** 2 / E2[k] if E2[k] > 0 else 0.
        rows.append(row)
        print(f"M={M} L={L}: presel {S['presel']:.4g} fb | L3 {S['L3']:.3g} | L4 {S['L4']:.3g} (N={N['L4']}) | endpoints L4 {nep}: chi2>3 {nep3} ({nep3/max(nep,1):.1%}), chi2>9 {nep9}, R<0.5m {nepR} | events chi2>3 {S['C3']:.3g} fb (N={N['C3']}, Neff={row['Neff_C3']:.1f}), both-endpoints {S['C3x2']:.3g} (N={N['C3x2']})", flush=True)
with open(OUT, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
