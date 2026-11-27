"""Map-based recast of the ATLAS disappearing-track search (ATLAS-SUSY-2018-19 = arXiv:2201.02472, EWK signal region,
model-independent selection) for decaying charged quirks simulated with the patched QuirkPheno transport.

The selection and the public efficiency maps follow the recast of Ramos and Lessa (llprecasting/recastingCodes,
DisappearingTracks/ATLAS-SUSY-2018-19, which uses the ATLAS HEPData efficiency maps and the model-independent upper limit
sigma_vis < 0.037 fb): MET trigger efficiency map eff(MET, pT(jet1)); MET > 200 GeV, pT(jet1) > 100 GeV, min dphi(jets, MET) > 1;
tracklet efficiency map eff(eta_LLP, R_decay); smeared tracklet pT > 60 GeV; 0.1 < |eta| < 1.9; dR(LLP, jet) > 0.4;
event efficiency = trigger x (1 - prod_i (1 - eff_i)) over the decayed charged LLPs.  Differences, all deliberate:
  * jets/MET are generator-level (anti-kT 0.4 from visible particles, jetmet_gen.py) with a simple parametrised smearing instead of Delphes;
  * the charged quirk decays inside the Geant4 simulation (R_decay, eta, pT, from the simulated decay), no lifetime reweighting;
  * the string changes the curvature of the track: the transverse bending angle accumulated inside r < 125 mm is compared with the
    magnetic expectation (B of the simulation, 3.8 T) and the string-induced part is converted into an equivalent shift of q/pT at
    the ATLAS field (2 T) that is added to the q/pT smearing of the tracklet.  delta_phi_s = 0 gives the ideal-helix (wino-like) efficiency.
Samples: recoil-enriched Pythia events (rec_gen.sh: generator-level pair pT > 180 GeV, |eta| < 2.5 as in the main scan), re-simulated with the same
patched transport (out3/triplet_{dy,w}_M<M>/L<Lambda>_c*.root); the normalisation uses the number of generated events before the pair-pT filter.
Only the triplet is analysed: the doublet (c*tau ~ 7-10 mm) decays inside the first pixel layers, where the ATLAS tracklet map is zero.
usage: atlas_dt_recast.py <ms dir> <ATLAS_data dir> <kfactor_check.csv> <out.csv> [smear=first|intended]
smear: 'first' reproduces the pT-bin selection of the recast code (first bin that matches, i.e. alpha=1.86, sigma=20.94 /TeV for every pT > 10 GeV),
'intended' uses the highest pT bin below the track pT."""
import sys, re, csv, glob, numpy as np, uproot
base, atl, kf, out = sys.argv[1:5]
SMEAR = sys.argv[5] if len(sys.argv) > 5 else 'first'
import os
HITFIT = os.environ.get('HITFIT', '0') == '1'; HITVETO = os.environ.get('HITVETO', '0') == '1'
MASSES = [300, 500, 750, 1000, 1500, 2000]; LAMS = [30, 100, 300, 1000]
KW2 = {'doublet': 1., 'triplet': 2.}; SPECIES = {'doublet': 1., 'triplet': 2.}
BSIM, BATL = 3.8, 2.0
rng = np.random.default_rng(2024)

# ---------------- ATLAS maps
f = uproot.open(f'{atl}/DisappearingTrack2018-EfficiencyMaps.root')
def hist2(n):
    h = f[n]; return h.values(), h.axes[0].edges(), h.axes[1].edges()
TRIG, TX, TY = hist2('eff_trigger_average'); TRIG = np.nan_to_num(TRIG)
TRK, KX, KY = hist2('h_effmap_average_EWK'); TRK = np.nan_to_num(TRK)
def look(v, ex, ey, x, y):
    i = np.searchsorted(ex, x, side='right') - 1; j = np.searchsorted(ey, y, side='right') - 1
    if i < 0 or j < 0 or i >= v.shape[0] or j >= v.shape[1]: return 0.
    return float(v[i, j])
# ---------------- q/pT smearing (crystal-ball-like core + exponential tails; alpha, sigma in 1/TeV)
PAIRS = [(10., 1.86, 20.94), (15., 1.86, 19.54), (20., 1.86, 18.33), (25., 1.86, 17.01), (35., 1.82, 15.42), (45., 1.66, 14.49), (60., 1.54, 13.90), (100., 1.64, 14.03)]
def make_sampler(alpha, sigma):
    x = np.linspace(-800., 800., 32001); z = x / sigma
    pdf = np.where(z < -alpha, np.exp(alpha * (z + alpha / 2)), np.where(z > alpha, np.exp(-alpha * (z - alpha / 2)), np.exp(-z ** 2 / 2)))
    cdf = np.cumsum(pdf); cdf /= cdf[-1]
    return lambda u: np.interp(u, cdf, x)
SAMP = [(p, make_sampler(a, s)) for p, a, s in PAIRS]
def pick(pt):
    if pt < 10.: return None
    if SMEAR == 'first': return SAMP[0][1]
    s = None
    for p, fn in SAMP:
        if pt > p: s = fn
    return s
# ---------------- jets / MET smearing (generic; stated in the paper)
def jet_sigma(pt): return np.sqrt(0.5 ** 2 / np.maximum(pt, 1.) + 0.03 ** 2)
SOFT = 8.0   # GeV per MET component from unclustered energy / pile-up-free soft term
def dphi(a, b): return abs((a - b + np.pi) % (2 * np.pi) - np.pi)
def event_objects(jm, k):
    jets = jm['jets'][k]; jets = jets[jets[:, 0] > 0]
    pt, eta, ph = jets[:, 0].copy(), jets[:, 1], jets[:, 2]
    mx, my = jm['metx'][k], jm['mety'][k]
    fct = 1. + rng.normal(0., jet_sigma(pt))
    # replace the true jet vectors by the smeared ones in the MET
    dpt = pt * (fct - 1.); mx -= np.sum(dpt * np.cos(ph)); my -= np.sum(dpt * np.sin(ph))
    mx += rng.normal(0., SOFT); my += rng.normal(0., SOFT)
    pt = pt * fct; sel = (np.abs(eta) < 2.8) & (pt > 20.)
    o = np.argsort(-pt[sel])
    return pt[sel][o], eta[sel][o], ph[sel][o], float(np.hypot(mx, my)), float(np.arctan2(my, mx))
def presel(pt, eta, ph, met, mph):
    """returns (trigger eff, kinematic pass) following the recast code; jets are pT>20, |eta|<2.8"""
    if len(pt) == 0 or np.any(np.abs(eta) > 2.4): return 0., False
    trig = 1. if met > 500. else look(TRIG, TX, TY, met, min(499., pt[0]))
    kin = met >= 200. and pt[0] >= 100.
    if kin:
        d = [dphi(mph, ph[i]) for i in range(min(4, len(pt))) if pt[i] >= 50.]
        kin = (min(d) if d else 99999.) > 1.0
    return trig, kin
# ---------------- K-factors (public NLO+NLL / own LO), log-linear in M
kr = np.genfromtxt(kf, delimiter=',', names=True)
def K(M, ch): return float(np.interp(np.log(M), np.log(kr['M_GeV']), kr['K_DY_own' if ch == 'dy' else 'K_W_own']))
def sigma_pb(fn): return float(re.search(r'sigma = ([\d.eE+-]+) pb', open(fn).read()).group(1))
# ---------------- main loop
def quirks(a, i, Lm=125.):
    """decayed charged quirks of event i: dict(pT, eta, phi, R, dphis, Lm, sq)"""
    res = []
    for j, s in enumerate(a['dec_side'][i]):
        s = int(s); pdg = int(a['q_pdg'][i][s])
        if abs(pdg) != 17: continue
        ppx, ppy, ppz = a['dec_ppx'][i][j] / 1e3, a['dec_ppy'][i][j] / 1e3, a['dec_ppz'][i][j] / 1e3
        pt = np.hypot(ppx, ppy); eta = np.arcsinh(ppz / max(pt, 1e-9)); R = np.hypot(a['dec_x'][i][j], a['dec_y'][i][j])
        sel = (a['h_q'][i] == s) & np.isin(a['h_det'][i], [2, 4, 5, 6, 7, 8]); dphis = 0.; L = 0.1; chi2ok = 1.
        if sel.sum() >= 2:
            o = np.argsort(a['h_t0'][i][sel])
            x = ((a['h_x0'][i] + a['h_x1'][i]) / 2)[sel][o]; y = ((a['h_y0'][i] + a['h_y1'][i]) / 2)[sel][o]
            r = np.hypot(x - a['q_x0'][i][s], y - a['q_y0'][i][s]); m = np.where(r < Lm)[0]
            if len(m) >= 2:
                jj = m[-1]; px = a['h_px'][i][sel][o][jj]; py = a['h_py'][i][sel][o][jj]
                p0x, p0y = a['q_px0'][i][s], a['q_py0'][i][s]; pt0 = np.hypot(p0x, p0y) / 1e3; Lm_ = r[jj] / 1e3
                d = (np.arctan2(py, px) - np.arctan2(p0y, p0x) + np.pi) % (2 * np.pi) - np.pi
                sq = 1. if pdg > 0 else -1.
                dphis = d - sq * 0.29979 * BSIM * Lm_ / pt0; L = Lm_
                if HITFIT and len(m) >= 2:
                    # hit-level alternative: circle through the production point (beam spot) and the true hit positions at r < Lm, signed
                    # curvature relative to the initial direction; the string-induced part is kappa_fit - kappa_B, same conversion afterwards
                    xs_, ys_ = x[m] - a['q_x0'][i][s], y[m] - a['q_y0'][i][s]
                    A = np.c_[xs_, ys_]; bvec = 0.5 * (xs_ ** 2 + ys_ ** 2)
                    ctr = np.linalg.lstsq(A, bvec, rcond=None)[0]; R2 = ctr @ ctr
                    ph0 = np.arctan2(p0y, p0x)
                    kap = (-ctr[0] * np.sin(ph0) + ctr[1] * np.cos(ph0)) / R2 * 1e3 if R2 > 0 else 0.    # 1/m, CCW positive
                    # chi2-like non-helicity: rms residual of hits to the fitted circle (mm)
                    dphis = (kap - sq * 0.29979 * BSIM / pt0) * Lm_
                    # deterministic non-helicity veto: string-induced hit residuals (no hit noise) against sigma = 10 um; ATLAS requires
                    # chi2 probability > 0.1 (chi2 > 2.71 for one degree of freedom, n_hits = 3 here) and no significantly deviating hit
                    if HITVETO and len(m) >= 3:
                        rsd = np.hypot(xs_ - ctr[0], ys_ - ctr[1]) - np.sqrt(R2)
                        chi2ok = float(np.sum((rsd / 0.010) ** 2) < 2.71)
        res.append(dict(chi2ok=chi2ok, pt=pt, pt0=np.hypot(a['q_px0'][i][s], a['q_py0'][i][s]) / 1e3, eta=eta, phi=np.arctan2(ppy, ppx), R=R, dphis=dphis, L=L, sq=(1. if pdg > 0 else -1.)))
    return res
def track_eff(q, jets_eta, jets_phi, string):
    eff = look(TRK, KX, KY, q['eta'], q['R'])
    if eff <= 0.: return 0.
    if string: eff *= q['chi2ok']
    qpt = q['sq'] * 1e3 / q['pt0']
    if string and q['L'] > 0: qpt += 1e3 * q['dphis'] / (0.29979 * BATL * q['L'])
    smp = pick(q['pt0'])
    if smp is None: return 0.
    ptm = 1e3 / abs(qpt + smp(q['u']))
    if ptm <= 60. or not (0.1 < abs(q['eta']) < 1.9): return 0.
    for e, p in zip(jets_eta, jets_phi):
        if np.hypot(q['eta'] - e, dphi(q['phi'], p)) < 0.4: return 0.
    return eff
params = {}
for l in open(f'{base}/ms_params.txt'):
    m = re.match(r'M=(\d+) (\w+) dm=([\d.eE+-]+) GeV ctau=([\d.]+) mm', l); params[(int(m[1]), m[2])] = (float(m[3]), float(m[4]))
MASSES = [300, 500, 750, 1000, 1500]; rep = 'triplet'
ESCS = {1, 2}
rows = []
for M in MASSES:
    sg, ng = {}, {}
    for p in ('dy', 'w'):
        t = open(f'{base}/M{M}/rec_{p}.log').read(); sg[p] = float(re.search(r'sigma = ([\d.eE+-]+) pb', t).group(1)) * 1e3; ng[p] = int(re.search(r'wrote (\d+) events', t).group(1))
    rw = np.load(f'{base}/M{M}/rec_dy_rw.npz'); jmz = {p: np.load(f'{base}/M{M}/jm3_{p}.npz') for p in ('dy', 'w')}
    for L in LAMS:
        keys = ('vis', 'vis_ideal', 'presel', 'pair_presel', 'sig'); acc = {k: 0. for k in keys}; accNLO = {k: 0. for k in keys}; err2 = {k: 0. for k in keys}
        for p in ('dy', 'w'):
            fn = glob.glob(f'{base}/out3/{rep}_{p}_M{M}/L{L}_c*.root')[0]
            a = uproot.open(fn)['events'].arrays(['genEvt', 'q_pdg', 'q_reason', 'q_x0', 'q_y0', 'q_px0', 'q_py0', 'dec_side', 'dec_x', 'dec_y', 'dec_ppx', 'dec_ppy', 'dec_ppz',
                                                 'h_q', 'h_det', 'h_x0', 'h_y0', 'h_x1', 'h_y1', 'h_px', 'h_py', 'h_t0'], library='np')
            ids = a['genEvt'].astype(int); scale = sg[p] / ng[p]
            if p == 'dy': wmap = dict(zip(rw['evnum'], rw[rep])); w = SPECIES[rep] * np.array([wmap[i] for i in ids])
            else: w = np.full(len(ids), KW2[rep] * SPECIES[rep])
            jm = jmz[p]; idx = {int(e): k for k, e in enumerate(jm['evnum'])}
            s = {k: 0. for k in keys}; s2 = {k: 0. for k in keys}; kk = K(M, p)
            for i, g in enumerate(ids):
                pt, eta, ph, met, mph = event_objects(jm, idx[g])
                trig, kin = presel(pt, eta, ph, met, mph)
                if trig <= 0. or not kin: continue
                rs = [int(v) for v in a['q_reason'][i]]; ndec = len(a['dec_side'][i])
                esc = all(v in ESCS or v == 4 for v in rs) and any(v in ESCS for v in rs)
                qs = quirks(a, i)
                for q in qs: q['u'] = rng.random()
                contrib = dict(presel=w[i] * trig, pair_presel=w[i] * trig * float(ndec > 0 and esc))
                for key, st in (('vis', True), ('vis_ideal', False)):
                    pr = 1.
                    for q in qs: pr *= 1. - track_eff(q, eta, ph, st)
                    contrib[key] = w[i] * trig * (1. - pr)
                for key, v in contrib.items(): s[key] += v; s2[key] += v * v
            for key in s: acc[key] += scale * s[key]; accNLO[key] += scale * s[key] * kk; err2[key] += (scale * kk) ** 2 * s2[key]
        dm, ct = params[(M, rep)]
        rows.append(dict(M_GeV=M, rep=rep, Lambda_eV=L, ctau_mm=ct, sigma_vis_LO_fb=acc['vis'], sigma_vis_NLO_fb=accNLO['vis'], sigma_vis_NLO_err_fb=err2['vis'] ** .5,
                         sigma_vis_idealhelix_NLO_fb=accNLO['vis_ideal'], sigma_vis_idealhelix_NLO_err_fb=err2['vis_ideal'] ** .5,
                         sigma_presel_NLO_fb=accNLO['presel'], sigma_pair_presel_NLO_fb=accNLO['pair_presel'], sigma_pair_presel_NLO_err_fb=err2['pair_presel'] ** .5))
        print(f"M={M} L={L}: sigma_vis NLO {accNLO['vis']:.4g} +- {err2['vis']**.5:.2g} (ideal helix {accNLO['vis_ideal']:.4g} +- {err2['vis_ideal']**.5:.2g}) fb, presel {accNLO['presel']:.4g}, neutral pair + presel {accNLO['pair_presel']:.4g}", flush=True)
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
