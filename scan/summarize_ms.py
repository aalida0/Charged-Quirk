"""Combine the mass-scan outputs into yields: M x {doublet, triplet} x Lambda, DY (re-weighted) + W*-mediated production.
usage: summarize_ms.py <ms dir containing M<mass>/, csv/ and out/> <out.csv>
sigma_eff = sigma_gen(Pythia) * [sum_i w_i flag_i] / (max event id covered + 1);  DY: w = LO coupling re-weighting (reweight_dy.py);
W*: w = kappa_W^2 (1 doublet, 2 triplet: |<1,1|T+|1,0>|^2 = 2); both: x2 charged Dirac species for the triplet (Psi+ and Psi-).  |eta| < 2.5 generator filter; no k-factor."""
import sys, re, csv, numpy as np
base, out = sys.argv[1], sys.argv[2]
MASSES = [300, 500, 750, 1000, 1500, 2000]; LAMS = [30, 100, 300, 1000]
KW2 = {'doublet': 1., 'triplet': 2.}      # W coupling^2 relative to the doublet
SPECIES = {'doublet': 1., 'triplet': 2.}  # charged Dirac species: doublet (Y=1/2) Psi+; triplet (Y=0) Psi+ and Psi- (a Dirac triplet = 2 Majorana winos)
def sigma_pb(f): return float(re.search(r'sigma = ([\d.eE+-]+) pb', open(f).read()).group(1))
params = {}
for l in open(f'{base}/ms_params.txt'):
    m = re.match(r'M=(\d+) (\w+) dm=([\d.eE+-]+) GeV ctau=([\d.]+) mm', l); params[(int(m[1]), m[2])] = (float(m[3]), float(m[4]))
rows = []
for M in MASSES:
    sg = {p: sigma_pb(f'{base}/M{M}/{p}.log') * 1e3 for p in ('dy', 'w')}      # fb
    rw = np.load(f'{base}/M{M}/dy_rw.npz')
    for rep in ('doublet', 'triplet'):
        dm, ct = params[(M, rep)]
        for L in LAMS:
            tot = dict(term=0., neut=0., metv=0., tout=0.); sig = {}
            for p in ('dy', 'w'):
                ev = np.load(f'{base}/csv/{rep}_{p}_M{M}_ev.npz')[f'L{L}']
                ids = ev[:, 0].astype(int); cover = ids.max() + 1
                if p == 'dy':
                    wmap = dict(zip(rw['evnum'], rw[rep])); w = SPECIES[rep] * np.array([wmap[i] for i in ids])
                else: w = np.full(len(ids), KW2[rep] * SPECIES[rep])
                scale = sg[p] / cover
                sig[p] = scale * w.sum()
                for k, c in (('term', 1), ('neut', 2), ('metv', 3), ('tout', 4)): tot[k] += scale * (w * ev[:, c]).sum()
            h = []
            for p in ('dy', 'w'):
                hx = np.genfromtxt(f'{base}/csv/helix_{rep}_{p}_M{M}.csv', delimiter=',', names=True)
                hx = hx[np.isclose(hx['Lambda_eV'], L)][0]
                h.append((float(hx['n_ge4hits']), float(hx['f_helical'])))
            nh = sum(x[0] for x in h); fh = sum(x[0] * x[1] for x in h if x[0] > 0) / nh if nh > 0 else np.nan
            rows.append(dict(M_GeV=M, rep=rep, Lambda_eV=L, dm_MeV=dm * 1e3, ctau_mm=ct, sigma_DY_fb=sig['dy'], sigma_Wstar_fb=sig['w'],
                sigma_tot_fb=sig['dy'] + sig['w'], sigma_termtrack_fb=tot['term'], sigma_neutralpair_escape_fb=tot['neut'],
                sigma_neutralpair_MET250_fb=tot['metv'], sigma_timeout_fb=tot['tout'], sigma_DTlike_fb=tot['term'] * (fh if fh == fh else 0.), n_quirks_ge4hits=nh, f_helical_ge4hits=fh))
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows:
    if r['Lambda_eV'] == 100: print(f"{r['rep']:7s} M={r['M_GeV']:5d} dm={r['dm_MeV']:.0f} ct={r['ctau_mm']:.1f}mm sigma DY/W*={r['sigma_DY_fb']:.3g}/{r['sigma_Wstar_fb']:.3g} fb term={r['sigma_termtrack_fb']:.3g} fb nMET={r['sigma_neutralpair_MET250_fb']:.3g} fb")
