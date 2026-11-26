#!/usr/bin/env python3
"""Validate a quirksim run with the decay step on (gun pair, one fixed (m, Delta, c tau)).
usage: validate_decay.py out.root <M_GeV> <Delta_GeV> <ctau_mm>
Checks: (1) proper decay time ~ Exp(c tau); (2) 4-momentum conservation per decay; (3) rest-frame product spectra vs
analytic phase space (pi: fixed |p*|; e, mu: E p E_nu^2 weighting of a constant-matrix-element 3-body decay)."""
import sys
import numpy as np
import uproot
from scipy import stats
from scipy.integrate import quad

f, M, D, ctau = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
M0 = M - D
mpi, me, mmu = 0.13957, 0.000511, 0.10566
a = uproot.open(f)['events'].arrays(['dec_side', 'dec_t', 'dec_pe', 'dec_ppx', 'dec_ppy', 'dec_ppz', 'dec_de', 'dec_dpx',
                                     'dec_dpy', 'dec_dpz', 'dp_rec', 'dp_pdg', 'dp_px', 'dp_py', 'dp_pz', 'dp_e', 'q_dec'], library='np')
# ntuple units: mm, ns, MeV -> convert momenta/energies to GeV
for k in a:
    if k.startswith(('dec_p', 'dec_d', 'dp_p', 'dp_e')) and k not in ('dp_pdg', 'dp_rec'):
        a[k] = [np.asarray(x) * 1e-3 for x in a[k]]
n_ev = len(a['dec_side'])
n_dec = sum(len(x) for x in a['dec_side'])
print(f'events {n_ev}, decays {n_dec} (2 per event expected: {2*n_ev}); quirks flagged decayed: {int(sum(np.sum(x) for x in a["q_dec"]))}')

# (1) lifetime: c * t_lab / gamma, gamma = E/M at decay (the quirk moves almost freely over millimetres)
c_mm_ns = 299.792458
tau_mm = np.concatenate([c_mm_ns * np.asarray(t) * M / np.asarray(pe) for t, pe in zip(a['dec_t'], a['dec_pe'])])
print(f'(1) proper c*tau: mean {tau_mm.mean():.3f} +- {tau_mm.std()/np.sqrt(len(tau_mm)):.3f} mm (input {ctau}); '
      f'KS vs Exp({ctau}): p = {stats.kstest(tau_mm, "expon", args=(0, ctau)).pvalue:.3f}')

# (2) conservation and (3) rest-frame spectra
worst = 0.; pistar, ee, emu = [], [], []
for ie in range(n_ev):
    rec = np.asarray(a['dp_rec'][ie]).astype(int)
    for ir in range(len(a['dec_side'][ie])):
        P = np.array([a['dec_pe'][ie][ir], a['dec_ppx'][ie][ir], a['dec_ppy'][ie][ir], a['dec_ppz'][ie][ir]])
        Dd = np.array([a['dec_de'][ie][ir], a['dec_dpx'][ie][ir], a['dec_dpy'][ie][ir], a['dec_dpz'][ie][ir]])
        sel = np.nonzero(rec == ir)[0]
        prods = [(int(a['dp_pdg'][ie][k]), np.array([a['dp_e'][ie][k], a['dp_px'][ie][k], a['dp_py'][ie][k], a['dp_pz'][ie][k]])) for k in sel]
        tot = Dd + sum((p for _, p in prods), np.zeros(4))
        worst = max(worst, np.abs(tot - P).max())
        # boost to the parent rest frame
        b = P[1:] / P[0]; g = P[0] / M; b2 = b @ b
        def rest(p):
            bp = b @ p[1:]
            return np.concatenate([[g * (p[0] - bp)], p[1:] + ((g - 1) * bp / b2 - g * p[0]) * b])
        pdgs = [p for p, _ in prods]
        if pdgs and len(prods) == 1 and abs(pdgs[0]) == 211:
            pistar.append(np.linalg.norm(rest(prods[0][1])[1:]))
        for pdg, p in prods:
            if abs(pdg) == 11: ee.append(rest(p)[0])
            if abs(pdg) == 13: emu.append(rest(p)[0])
print(f'(2) max |sum(final) - parent| over components: {worst*1e6:.3g} keV (GeV-scale momenta; float precision of the ntuple)')
lam = (M**2 - (M0 + mpi)**2) * (M**2 - (M0 - mpi)**2)
print(f'(3) pi- alone: rest-frame |p*| mean {np.mean(pistar)*1e3:.4f} MeV, spread {np.std(pistar)*1e3:.4f} MeV; analytic {np.sqrt(lam)/(2*M)*1e3:.4f} MeV  (N={len(pistar)})')
def mean_e(m, Emax):
    w = lambda E: E * np.sqrt(max(E*E - m*m, 0)) * (Emax - E)**2
    return quad(lambda E: E * w(E), m, Emax)[0] / quad(w, m, Emax)[0]
for name, arr, m in (('e', ee, me), ('mu', emu, mmu)):
    if len(arr) > 20:
        Emax = (M**2 - M0**2 + m*m) / (2*M)
        print(f'    {name}: N={len(arr)}, mean rest-frame E {np.mean(arr)*1e3:.2f} +- {np.std(arr)/np.sqrt(len(arr))*1e3:.2f} MeV; '
              f'analytic E p (Emax-E)^2 mean {mean_e(m, Emax)*1e3:.2f} MeV (Emax {Emax*1e3:.2f})')
