"""C6 check of the transport proper-time bookkeeping.  Every charged quirk eventually decays in these samples (no censoring), so for a
constant-lifetime exponential law the mean proper decay time must equal ctau0.  Proper time at the decay is estimated as
t_lab / <gamma>, <gamma> = mean of the start and decay-point gamma (exact for negligible string loss, Lambda=30 eV; an approximation
that can bias the result by a few % when the string slows the quirk).  Also a one-sample KS test against Exp(ctau0).
usage: proper_time_check.py <dir with L<Lambda>_c<ctau>.root> <M_GeV>"""
import sys, glob, re, numpy as np, uproot
d, M = sys.argv[1], float(sys.argv[2]) * 1e3
c = 299.792458  # mm/ns
for f in sorted(glob.glob(f'{d}/L*_c*.root'), key=lambda x: float(re.search(r'/L([\d.]+)_', x).group(1))):
    L = float(re.search(r'/L([\d.]+)_', f).group(1)); ct0 = float(re.search(r'_c([\d.]+)\.root', f).group(1))
    a = uproot.open(f)['events'].arrays(['q_pdg', 'q_e0', 'dec_side', 'dec_t', 'dec_pe'], library='np')
    tau = []
    for i in range(len(a['q_pdg'])):
        for j, s in enumerate(a['dec_side'][i]):
            s = int(s)
            if abs(int(a['q_pdg'][i][s])) != 17: continue
            tau.append(a['dec_t'][i][j] * c / (0.5 * (a['q_e0'][i][s] + a['dec_pe'][i][j]) / M))
    tau = np.sort(tau); n = len(tau); F = 1 - np.exp(-tau / ct0)
    D = max(np.max(np.arange(1, n + 1) / n - F), np.max(F - np.arange(0, n) / n))
    print(f'L={L:6.0f} eV ctau0={ct0:7.3f} mm N={n} <c*tau>={tau.mean():.2f}+-{tau.std()/np.sqrt(n):.2f} mm ratio={tau.mean()/ct0:.3f} KS D={D:.4f} (5% crit {1.36/np.sqrt(n):.4f})')
