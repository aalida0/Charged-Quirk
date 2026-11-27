"""Generator-level recast of the ATLAS monojet search (arXiv:2102.10874, 139/fb, 13 TeV) for the quirk pair, all production channels,
against the published model-independent limits on sigma x A x eps of the inclusive regions IM0-IM10 (observed, in fb).
Selection (generator-level jets/MET of jetmet_gen.py, quirks invisible, same parametrised jet/MET smearing as atlas_dt_recast.py):
 MET > cut; leading jet pT > 150 GeV, |eta| < 2.4; at most 4 jets with pT > 30 GeV, |eta| < 2.8; Delta phi(jet_i, MET) > 0.4 (0.6 for MET <= 250) for those jets;
 lepton/photon veto trivially passed; trigger taken as fully efficient above 200 GeV.  The detector efficiency eps is NOT applied (eps = 1),
 so the fiducial cross section is an upper bound on sigma x A x eps; "R = sigma_fid / sigma_lim" < 1 means not constrained even for eps = 1.
Channels: Drell-Yan Q+Q- and W*-mediated Q+Q0 (as in the other recasts), and for the doublet the Z-mediated neutral pair Q0 Q0bar
(LO re-weighting of the same DY events, reweight_dy_n0.py), triplet: Q0 does not couple to the Z.  NLO+NLL factors: K_DY for Q+Q- and Q0Q0bar, K_W for the W* channel.
Recoil-enriched samples (generator pair pT > 180 GeV): valid for MET >= 250 GeV (events with pair pT < 180 GeV cannot reach MET > 250 GeV except through tails of the smearing).
usage: monojet_recast.py <ms dir> <kfactor_check.csv> <out.csv>"""
import sys, re, csv, numpy as np
base, kf, out = sys.argv[1:4]
rng = np.random.default_rng(7)
MASSES = [300, 500, 750, 1000, 1500]
KW2 = {'doublet': 1., 'triplet': 2.}; SPECIES = {'doublet': 1., 'triplet': 2.}
CUTS = [200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000]
LIM = {200: 736., 250: 296., 300: 150., 350: 86., 400: 52., 500: 21., 600: 10., 700: 4.1, 800: 2.1, 900: 1.2, 1000: 1.3}   # fb, observed, arXiv:2102.10874
SOFT = 8.0
def jet_sigma(pt): return np.sqrt(0.5 ** 2 / np.maximum(pt, 1.) + 0.03 ** 2)
def dphi(a, b): return abs((a - b + np.pi) % (2 * np.pi) - np.pi)
kr = np.genfromtxt(kf, delimiter=',', names=True)
def K(M, ch): return float(np.interp(np.log(M), np.log(kr['M_GeV']), kr['K_DY_own' if ch == 'dy' else 'K_W_own']))
def objects(jm, k):
    jets = jm['jets'][k]; jets = jets[jets[:, 0] > 0]; pt, eta, ph = jets[:, 0].copy(), jets[:, 1], jets[:, 2]
    mx, my = jm['metx'][k], jm['mety'][k]; fct = 1. + rng.normal(0., jet_sigma(pt)); dpt = pt * (fct - 1.)
    mx -= np.sum(dpt * np.cos(ph)); my -= np.sum(dpt * np.sin(ph)); mx += rng.normal(0., SOFT); my += rng.normal(0., SOFT)
    pt = pt * fct; sel = (np.abs(eta) < 2.8) & (pt > 30.); o = np.argsort(-pt[sel])
    return pt[sel][o], eta[sel][o], ph[sel][o], float(np.hypot(mx, my)), float(np.arctan2(my, mx))
def passes(pt, eta, ph, met, mph, cut):
    if met < cut or len(pt) == 0 or pt[0] < 150. or abs(eta[0]) > 2.4 or len(pt) > 4: return False
    dmin = 0.6 if met <= 250. else 0.4
    return all(dphi(mph, p) > dmin for p in ph)
def sigma_fb(fn): return float(re.search(r'sigma = ([\d.eE+-]+) pb', open(fn).read()).group(1)) * 1e3
rows = []
for M in MASSES:
    sg, ng = {}, {}
    for p in ('dy', 'w'):
        t = open(f'{base}/M{M}/rec_{p}.log').read(); sg[p] = float(re.search(r'sigma = ([\d.eE+-]+) pb', t).group(1)) * 1e3; ng[p] = int(re.search(r'wrote (\d+) events', t).group(1))
    rw = np.load(f'{base}/M{M}/rec_dy_rw.npz'); rn = np.load(f'{base}/M{M}/rec_dy_rw_n0.npz')
    jm = {p: np.load(f'{base}/M{M}/jm3_{p}.npz') for p in ('dy', 'w')}
    for rep in ('doublet', 'triplet'):
        chans = [('dy', 'dy', dict(zip(rw['evnum'], rw[rep])), SPECIES[rep]), ('w', 'w', None, KW2[rep] * SPECIES[rep])]
        if rep == 'doublet': chans.append(('dy0', 'dy', dict(zip(rn['evnum'], rn['doublet_n0'])), 1.))
        tot = {c: 0. for c in CUTS}; err = {c: 0. for c in CUTS}; per = {}
        for name, p, wmap, mult in chans:
            j = jm[p]; scale = sg[p] / ng[p] * K(M, p); acc = {c: 0. for c in CUTS}; acc2 = {c: 0. for c in CUTS}
            for k, e in enumerate(j['evnum']):
                w = mult * (wmap[int(e)] if wmap is not None else 1.)
                o = objects(j, k)
                for c in CUTS:
                    if passes(*o, c): acc[c] += w; acc2[c] += w * w
            for c in CUTS: tot[c] += scale * acc[c]; err[c] += (scale ** 2) * acc2[c]; per[(name, c)] = scale * acc[c]
        for c in CUTS:
            rows.append(dict(M_GeV=M, rep=rep, MET_cut=c, sigma_fid_NLO_fb=tot[c], err_fb=err[c] ** .5, sigma_lim_fb=LIM[c], R=tot[c] / LIM[c],
                             dy_fb=per[('dy', c)], w_fb=per[('w', c)], n0_fb=per.get(('dy0', c), 0.)))
        print(f"M={M} {rep}: " + ' '.join(f"{c}:{tot[c]:.3g}({tot[c]/LIM[c]:.2g})" for c in CUTS), flush=True)
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
