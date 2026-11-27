"""Compare the Pythia-quirk production normalisation (DY re-weighting; W* x kappa_W^2) with Pythia's own SUSY pure-wino
LO processes (wino_pythia_check.py output) and with public NLO+NLL wino numbers.
usage: kfactor_check.py wino_pythia.txt generator_summary.txt <dir with M<M>/dy_rw.npz> out.csv [xsec_chgchg.csv xsec_chgn1.csv]
With the two optional tables (LHC SUSY Cross Section WG wino NLO+NLL, 13 TeV, fb, as distributed with the ATLAS-SUSY-2018-19 recast
in the llprecasting/recastingCodes repository) the K-factors K_DY(M) = NLO+NLL / own LO and K_W(M) are added at every mass.
All numbers are per Majorana wino (a Dirac triplet = 2 winos), total cross section, no acceptance cut, LO, pSet 20."""
import sys, re, numpy as np, csv
wp, gs, base, out = sys.argv[1:5]
TAB = sys.argv[5:7]
W = {}
for l in open(wp):
    m = re.match(r'M=(\d+) (\d+) .*?: ([\d.eE+-]+) \+- ([\d.eE+-]+) fb', l)
    if m: W.setdefault(int(m[1]), {})[int(m[2])] = float(m[3])
G = {}
for l in open(gs):
    m = re.match(r'M=(\d+) (dy|w) quirkgen: .* sigma = ([\d.eE+-]+) pb', l)
    if m: G.setdefault(int(m[1]), {})[m[2]] = float(m[3]) * 1e3        # fb
def tab(f, M):
    d = np.genfromtxt(f, delimiter=',', names=True); return float(np.exp(np.interp(M, d['mass'], np.log(d['xsec']))))
rows = []
for M in sorted(W):
    w = np.load(f'{base}/M{M}/dy_rw.npz')['triplet'].mean()
    own_dy, own_w = G[M]['dy'] * w, G[M]['w'] * 2.0                    # per wino: DY*<w>, W* * kappa_W^2
    p_dy, p_w = W[M][1591], W[M][1571] + W[M][1572]
    rows.append(dict(M_GeV=M, own_DY_fb=own_dy, pythia_wino_chipchim_fb=p_dy, ratio_DY=own_dy / p_dy,
                     own_Wstar_fb=own_w, pythia_wino_chipmchi0_fb=p_w, ratio_Wstar=own_w / p_w, pythia_DY_over_W=p_dy / p_w))
    if TAB:
        t_dy, t_w = tab(TAB[0], M), tab(TAB[1], M)
        rows[-1].update(public_NLONLL_chipchim_fb=t_dy, public_NLONLL_chipmchi0_fb=t_w, K_DY_over_pythiaSUSY=t_dy / p_dy, K_W_over_pythiaSUSY=t_w / p_w, K_DY_own=t_dy / own_dy, K_W_own=t_w / own_w)
        print(f"   public NLO+NLL: chi+chi- {t_dy:.4g} fb, chi+-chi0 {t_w:.4g} fb;  K_DY = {t_dy / own_dy:.3f}, K_W = {t_w / own_w:.3f} (own LO); {t_dy / p_dy:.3f}, {t_w / p_w:.3f} (Pythia SUSY LO)")
    print(f"M={M}: DY own/Pythia-SUSY = {own_dy / p_dy:.3f}   W* own/Pythia-SUSY = {own_w / p_w:.3f}   DY/W*(Pythia SUSY) = {p_dy / p_w:.3f}")
# public NLO+NLL anchor: ATLAS 2106.01676, wino chi1+- chi2 0 at 850 GeV: 3.42 fb (charge-summed); Pythia LO wino at 850 by log-interpolation
a, b = W[750], W[1000]; sa, sb = a[1571] + a[1572], b[1571] + b[1572]
lo850 = np.exp(np.log(sa) + (np.log(sb) - np.log(sa)) * (850 - 750) / 250)
print(f"Pythia-SUSY LO chi+-chi0 (sum of charges) at 850 GeV (log-interp 750-1000): {lo850:.3f} fb;  K = 3.42/LO = {3.42 / lo850:.2f}")
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
