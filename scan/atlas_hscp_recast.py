"""Mass limits for a long-lived charged EW-multiplet quirk (Lambda <~ 100 eV, string negligible) from the public ATLAS HSCP search
(ATLAS-SUSY-2018-42, arXiv:2205.06013, 139/fb, HEPData ins2080541), using the *published* cross-section limit curves.

Method: for each ATLAS chargino lifetime (1, 4, 10, 30 ns, stable) the HEPData table gives the 95% CL observed and expected upper limit on
sigma(pp -> chi+chi-, chi+-chi0) [pb] vs chargino mass together with the signal cross section ATLAS assumed (wino, NLO+NLL).
A Dirac SU(2)-triplet quirk with the same mass and the same charged-particle kinematics has twice the wino production rate
(two Majorana winos per Dirac triplet; our LO generator reproduces Pythia's SUSY pure-wino rates to 0.89-0.96, see kfactor_check.py),
so sigma_quirk(M) = R * sigma_ATLAS-wino(M), R = 2 (triplet). The mass limit is where R*sigma_wino(M) crosses the ATLAS limit.
For lifetimes shorter than the 1 ns table, the ATLAS acceptance table (>=1 charged LLP with pT>120 GeV, |eta|<1.8, r_decay>500 mm) is
used to rescale the 1 ns limit: sigma_UL(tau) ~ sigma_UL(1 ns) * A(1 ns)/A(tau); this is only indicative (no efficiency map exists).
Same-kinematics assumption: the acceptance/efficiency of the quirk pair is taken equal to that of the chargino pair; spin/string
differences are not modelled (string forces ignored for Lambda <~ 100 eV).

usage: atlas_hscp_recast.py <HEPData-ins2080541-csv dir> out.csv"""
import sys, re, csv, numpy as np
D, OUT = sys.argv[1:3]
R = 2.0

def blocks(fn):
    """split a HEPData csv with several tables into list of (header, rows)"""
    out, cur, hdr = [], None, None
    for l in open(fn):
        l = l.strip()
        if l.startswith('#') or not l:
            if cur: out.append((hdr, cur)); cur, hdr = None, None
            continue
        if hdr is None: hdr, cur = l, []
        else: cur.append([float(x) if x not in ('--', '') else np.nan for x in re.split(r',(?=(?:[^"]*"[^"]*")*[^"]*$)', l)])
    if cur: out.append((hdr, cur))
    return out

def load(tag):
    b = blocks(f'{D}/CrossSectionLimit,Chargino{tag}.csv')
    e = np.array([r[:2] for r in b[0][1]]); o = np.array(b[1][1]); s = np.array(b[2][1])
    return dict(M=o[:, 0] * 1e3, obs=o[:, 1], exp=dict(zip(e[:, 0] * 1e3, e[:, 1])), sig=dict(zip(s[:, 0] * 1e3, s[:, 1])), Msig=s[:, 0] * 1e3, sigv=s[:, 1])

def cross(M, ratio):
    """largest mass where ratio(M) = sigma_sig / UL >= 1 (log-interpolated crossing); ratio decreasing with M"""
    above = ratio >= 1
    if above.all(): return M[-1], '>'
    if not above.any(): return M[0], '<'
    i = np.where(above)[0].max()
    if i == len(M) - 1: return M[-1], '>'
    # first crossing from above after last point that is above
    l0, l1 = np.log(ratio[i]), np.log(ratio[i + 1])
    return M[i] + (M[i + 1] - M[i]) * l0 / (l0 - l1), '='

rows = []
for tag, tau in [('1ns', 1.0), ('4ns', 4.0), ('10ns', 10.0), ('30ns', 30.0), ('Stable', np.inf)]:
    t = load(tag)
    M = t['M']; sig = np.array([t['sig'][m] for m in M]) if all(m in t['sig'] for m in M) else np.exp(np.interp(M, t['Msig'], np.log(t['sigv'])))
    obs = t['obs']; exp = np.array([t['exp'].get(m, np.nan) for m in M])
    for name, ul in [('observed', obs), ('expected', exp)]:
        ok = np.isfinite(ul)
        ratio_w = sig[ok] / ul[ok]                      # ATLAS wino
        mw, fw = cross(M[ok], ratio_w)
        mq, fq = cross(M[ok], R * ratio_w)              # Dirac triplet quirk
        rows.append(dict(tau_ns=tau, kind=name, wino_limit_GeV=round(mw), wino_flag=fw, triplet_quirk_limit_GeV=round(mq), quirk_flag=fq,
                         mass_range=f'{M[ok].min():.0f}-{M[ok].max():.0f}'))
        print(f"tau={tau} ns {name:8s}: wino {fw}{mw:.0f}  Dirac-triplet quirk (R={R:g}) {fq}{mq:.0f} GeV   (table covers {M[ok].min():.0f}-{M[ok].max():.0f})")

# below 1 ns: acceptance rescaling indicator
acc = blocks(f'{D}/Acceptance,Chargino.csv')
def accur(i):
    a = np.array(acc[i][1]); return a[np.isfinite(a[:, 1])]
a02, a1 = accur(0), accur(1)
print('\nbenchmark charged-state decay length (triplet, Delta m ~ 165 MeV): c*tau ~ 55 mm = 0.18 ns')
t1 = load('1ns')
for M in (300, 400, 500, 600, 700, 800):
    A02 = np.interp(M, a02[:, 0], a02[:, 1]); A1 = np.interp(M, a1[:, 0], a1[:, 1])
    ul_eff = np.interp(M, t1['M'], t1['obs']) * A1 / max(A02, 1e-4)
    sq = R * np.exp(np.interp(M, t1['Msig'], np.log(t1['sigv'])))
    print(f"  M={M}: A(0.2ns)={A02:.4f} (ATLAS table, rounded) A(1ns)={A1:.3f}; indicative sigma_UL(0.2 ns) >~ {ul_eff:.3g} pb vs sigma_quirk {sq:.3g} pb  -> {'excluded' if sq > ul_eff else 'not excluded'}")
with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
