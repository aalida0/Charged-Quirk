"""Calibrated recasts of published limits (flat-upper-limit approximation).
usage: recast_limits.py <massscan_summary.csv> <hscp_recast.csv> <generator_summary.txt> <out.csv>
Inputs (verified from the abstracts / text of the papers, see docs):
  ATLAS 2201.02472: pure-wino charginos excluded up to 660 GeV at tau = 0.2 ns (136 fb^-1) -- our triplet has c*tau = 54-59 mm = 0.18-0.20 ns
  CMS 2004.05153 (EXO-19-010): wino excluded below 474 GeV at tau = 0.2 ns (884 GeV at 3 ns)
  CMS 2410.09164 (EXO-18-002): Drell-Yan tau' (Q = 1e, 4th-generation lepton) excluded below 1.02 TeV (ionization method), 101 fb^-1
Method.  Everything is done in our own LO cross sections, so the k-factor cancels to first order.  The upper limit on sigma x A x eff is
assumed independent of the mass above the calibration mass (background-free plateau); this is the main approximation.
  Disappearing tracks: the Dirac triplet is 2 Majorana winos, sigma_wino-eq(M) = sigma_triplet(M)/2 (LO, DY + W*, c*tau-independent).  A model point
  is excluded if sigma_triplet(M) * eps_rel(Lambda, M) >= sigma_wino-eq(M_cal), with eps_rel = [DT-like rate](Lambda)/[DT-like rate](30 eV)
  (terminating-track yield x helical fraction, Lambda-dependence of the signal efficiency; ~1 for Lambda <= 300 eV).
  HSCP: tau' plateau UL = sigma_LO(Pythia DY tau', 4th-gen lepton, filter-corrected) at 1.02 TeV; compare with the surviving cross section."""
import sys, re, numpy as np, pandas as pd
from scipy.interpolate import CubicSpline
S = pd.read_csv(sys.argv[1]); H = pd.read_csv(sys.argv[2]); gen = open(sys.argv[3]).read(); out = sys.argv[4]
Ms = np.array(sorted(S.M_GeV.unique()), float)
def sig_tot(rep, M):
    s = S[(S.rep == rep) & (S.Lambda_eV == 100)].sort_values('M_GeV'); return np.exp(CubicSpline(Ms, np.log(s.sigma_tot_fb.values))(M))
def eps_rel(rep, L, M):
    r = []
    for m in (500., 750., 1000.):
        a = S[(S.rep == rep) & (S.M_GeV == m)]; r.append(float(a[a.Lambda_eV == L].sigma_DTlike_fb.iloc[0] / a[a.Lambda_eV == 30].sigma_DTlike_fb.iloc[0]))
    return float(np.interp(M, [500., 750., 1000.], r))
rows = []
grid = np.arange(250., 1800., 1.)
for lab, Mcal in (('ATLAS 2201.02472 (660 GeV)', 660.), ('CMS 2004.05153 (474 GeV)', 474.)):
    target = sig_tot('triplet', Mcal) / 2          # wino-equivalent LO at the calibration mass
    for L in (30, 100, 300, 1000):
        ex = [M for M in grid if sig_tot('triplet', M) * eps_rel('triplet', L, M) >= target]
        rows.append(dict(search=lab, rep='Dirac triplet', Lambda_eV=L, UL_fb=target, M_excl_GeV=max(ex) if ex else np.nan,
                         eps_rel_at_Mexcl=eps_rel('triplet', L, max(ex)) if ex else np.nan))
# HSCP: DY tau' plateau
def sig_dy_raw(M):
    m = [(int(a), float(b)) for a, b in re.findall(r'M=(\d+) dy quirkgen: wrote \d+ events, sigma = ([\d.eE+-]+) pb', gen)]
    mm, ss = zip(*sorted(m)); return float(np.exp(np.interp(M, mm, np.log(np.array(ss) * 1e3))))
ul = sig_dy_raw(1020.)
for rep in ('doublet', 'triplet'):
    for M in (300, 500, 750, 1000):
        s = H[(H.rep == rep) & (H.M_GeV == M) & (H.physical == 1) & (H.radius == 'tracker')].iloc[0]
        rows.append(dict(search='CMS 2410.09164 DY tau\' (1.02 TeV) plateau', rep=rep, Lambda_eV=np.nan, UL_fb=ul, M_excl_GeV=M,
                         eps_rel_at_Mexcl=s.sigma_HSCP_fb))   # eps_rel column holds sigma_HSCP surviving [fb] here
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows).round(4).to_string())
