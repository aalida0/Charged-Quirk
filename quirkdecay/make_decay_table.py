#!/usr/bin/env python3
"""Write the quirkdecay channel table and the proper decay length for a multiplet quirk.
usage: make_decay_table.py <dm_GeV> <kappa> <M_GeV> <out.txt> [--charge -1]
Prints c*tau [mm] (use with /quirk/decayCtau) and writes channels for the quirk (positive PDG code) of charge --charge.
Exclusive channels (e, mu, tau, pi, pi pi0); the caller must stay below Delta m ~ 0.8 GeV (hadronic continuum not modelled)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ctau'))
from ctau import channel_widths, hbar, c_cm

def main():
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    D, kappa, M, out = float(a[0]), float(a[1]), float(a[2]), a[3]
    charge = int(sys.argv[sys.argv.index('--charge') + 1]) if '--charge' in sys.argv else -1
    w = channel_widths(D, kappa, M0=M - D)
    tot = sum(w.values())
    s = '-' if charge < 0 else '+'
    # particle names for charge -1 (Geant4); lepton/pion names flip with the charge
    pi, pi0 = ('pi-' if charge < 0 else 'pi+'), 'pi0'
    lep = {'e': ('e-' if charge < 0 else 'e+', 'anti_nu_e' if charge < 0 else 'nu_e'),
           'mu': ('mu-' if charge < 0 else 'mu+', 'anti_nu_mu' if charge < 0 else 'nu_mu'),
           'tau': ('tau-' if charge < 0 else 'tau+', 'anti_nu_tau' if charge < 0 else 'nu_tau')}
    with open(out, 'w') as f:
        f.write(f'# charged quirk -> neutral quirk + X; Delta m = {D*1e3:.1f} MeV, kappa = {kappa}, M0 = {M-D:g} GeV\n')
        f.write('# weight = partial width / GeV (exclusive channels, static-source V-A current)\n')
        for k, v in w.items():
            if v <= 0: continue
            if k == 'pi': f.write(f'{v:.6e}  {pi}\n')
            elif k == 'pipi': f.write(f'{v:.6e}  {pi}  {pi0}\n')
            else: f.write(f'{v:.6e}  {lep[k][0]}  {lep[k][1]}\n')
    print(f'ctau_mm {hbar/tot*c_cm*10:.6g}')
    print('BR ' + ' '.join(f'{k}={v/tot:.4f}' for k, v in w.items()))

main()
