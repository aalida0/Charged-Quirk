"""Mass limits of a Dirac SU(2)-triplet quirk from atlas_dt_recast.py output (sigma_vis vs M) against the ATLAS-SUSY-2018-19 model-independent
limit sigma_vis < 0.037 fb, and validation of the same chain against the official ATLAS wino exclusion.
Validation: wino = half of the triplet (ideal-helix, i.e. no string) sigma_vis at the benchmark lifetime c*tau ~ 56 mm (tau = 0.19 ns, official
exclusion contour interpolated at that lifetime). The ratio  f = sigma_vis(M_official)/0.037 fb  quantifies how much our generator-level
proxy under- or over-estimates the visible cross section at the official limit. f < 1: the proxy is conservative (too small a visible cross section).
'calibrated' limits divide every sigma_vis by f (a single mass-independent factor, taken at the official limit mass), i.e. force the proxy to
reproduce the official wino limit at this lifetime; 'raw' limits are the uncalibrated proxy.
usage: atlas_dt_limits.py <atlas_dt_recast_*.csv> <official_exclusion.txt> <out.csv>"""
import sys, csv, numpy as np
fn, off, out = sys.argv[1:4]
UL = 0.037
d = np.genfromtxt(fn, delimiter=',', names=True)
o = np.array([[float(x) for x in l.split(';')] for l in open(off) if l.strip()])
tau = float(np.mean(d['ctau_mm'])) / 299.792458                # ns
M_off = float(np.interp(np.log(tau), np.log(o[:46, 1]), o[:46, 0]))   # lower (short-lifetime) branch of the contour
def sv(col, lam, M):
    m = d['Lambda_eV'] == lam; x, y = d['M_GeV'][m], d[col][m]
    y = np.maximum(y, 1e-12); return float(np.exp(np.interp(np.log(M), np.log(x), np.log(y))))
def lim(col, lam, f=1., fac=1.0, ecol=None):
    m = d['Lambda_eV'] == lam; x = d['M_GeV'][m]; y = fac * f * d[col][m]
    if ecol: y = y + 0                                    # (error handled by caller via shifted fac)
    if y.max() < UL: return x[0], '<'
    if y.min() >= UL: return x[-1], '>'
    i = np.where(y >= UL)[0].max()
    if i == len(x) - 1: return x[-1], '>'
    # use log-linear interpolation between the last point above and first point below the limit (y positive)
    y0, y1 = max(y[i], 1e-12), max(y[i + 1], 1e-12)
    return x[i] + (x[i + 1] - x[i]) * np.log(y0 / UL) / np.log(y0 / y1), '='
# ---- validation (wino = triplet/2, ideal helix, any Lambda)
lam0 = 30
w_at_off = 0.5 * sv('sigma_vis_idealhelix_NLO_fb', lam0, M_off)
f = w_at_off / UL
mw, fl = lim('sigma_vis_idealhelix_NLO_fb', lam0, fac=0.5)
print(f"validation: tau = {tau:.3f} ns; official wino lower-branch limit at this lifetime M = {M_off:.0f} GeV")
print(f"  proxy (ideal helix, wino = triplet/2): limit {fl}{mw:.0f} GeV; sigma_vis(M_off) = {w_at_off:.4f} fb vs 0.037 fb -> f = {f:.2f}  (f>1: proxy is optimistic at the official limit, f<1: conservative)")
rows = []
# N_IC: the tabulated sigma_vis is PER INFRACOLOR INDEX (Pythia samples carry no infracolor multiplicity); the Born rate of a Dirac multiplet in the
# fundamental of SU(N_IC) is N_IC times larger. The wino validation above does not involve N_IC.
for nic in (1, 2, 3):
  for lam in sorted(set(d['Lambda_eV'])):
    for name, col in [('string', 'sigma_vis_NLO_fb'), ('ideal_helix', 'sigma_vis_idealhelix_NLO_fb')]:
        mr, fr = lim(col, lam, fac=float(nic))
        mc, fc = lim(col, lam, f=1. / f if f > 0 else 1., fac=float(nic))
        rows.append(dict(N_IC=nic, Lambda_eV=lam, treatment=name, limit_raw_GeV=round(mr), flag_raw=fr, limit_calibrated_GeV=round(mc), flag_cal=fc, f=round(f, 3),
                         official_wino_limit_GeV=round(M_off), source=fn))
        print(f"N_IC={nic} Lambda={lam:5.0f} eV {name:12s}: raw {fr}{mr:.0f} GeV   calibrated to official wino {fc}{mc:.0f} GeV")
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
