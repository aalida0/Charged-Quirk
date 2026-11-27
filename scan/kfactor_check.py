"""Compare the Pythia-quirk production normalisation (DY re-weighting; W* x kappa_W^2) with Pythia's own SUSY pure-wino
LO processes (wino_pythia_check.py output) and with public NLO+NLL wino numbers.
usage: kfactor_check.py wino_pythia.txt generator_summary.txt <dir with M<M>/dy_rw.npz> out.csv
All numbers are per Majorana wino (a Dirac triplet = 2 winos), total cross section, no acceptance cut, LO, pSet 20."""
import sys, re, numpy as np, csv
wp, gs, base, out = sys.argv[1:5]
W = {}
for l in open(wp):
    m = re.match(r'M=(\d+) (\d+) .*?: ([\d.eE+-]+) \+- ([\d.eE+-]+) fb', l)
    if m: W.setdefault(int(m[1]), {})[int(m[2])] = float(m[3])
G = {}
for l in open(gs):
    m = re.match(r'M=(\d+) (dy|w) quirkgen: .* sigma = ([\d.eE+-]+) pb', l)
    if m: G.setdefault(int(m[1]), {})[m[2]] = float(m[3]) * 1e3        # fb
rows = []
for M in sorted(W):
    w = np.load(f'{base}/M{M}/dy_rw.npz')['triplet'].mean()
    own_dy, own_w = G[M]['dy'] * w, G[M]['w'] * 2.0                    # per wino: DY*<w>, W* * kappa_W^2
    p_dy, p_w = W[M][1591], W[M][1571] + W[M][1572]
    rows.append(dict(M_GeV=M, own_DY_fb=own_dy, pythia_wino_chipchim_fb=p_dy, ratio_DY=own_dy / p_dy,
                     own_Wstar_fb=own_w, pythia_wino_chipmchi0_fb=p_w, ratio_Wstar=own_w / p_w, pythia_DY_over_W=p_dy / p_w))
    print(f"M={M}: DY own/Pythia-SUSY = {own_dy / p_dy:.3f}   W* own/Pythia-SUSY = {own_w / p_w:.3f}   DY/W*(Pythia SUSY) = {p_dy / p_w:.3f}")
# public NLO+NLL anchor: ATLAS 2106.01676, wino chi1+- chi2 0 at 850 GeV: 3.42 fb (charge-summed); Pythia LO wino at 850 by log-interpolation
a, b = W[750], W[1000]; sa, sb = a[1571] + a[1572], b[1571] + b[1572]
lo850 = np.exp(np.log(sa) + (np.log(sb) - np.log(sa)) * (850 - 750) / 250)
print(f"Pythia-SUSY LO chi+-chi0 (sum of charges) at 850 GeV (log-interp 750-1000): {lo850:.3f} fb;  K = 3.42/LO = {3.42 / lo850:.2f}")
with open(out, 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
