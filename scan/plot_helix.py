"""f_helical vs Lambda from analyze_helix.py output.  usage: plot_helix.py helix_scan.csv out.png [helix_smeared.csv]  (solid: unsmeared RMS<0.1 mm; dashed: smeared-hit chi2 fit)"""
import sys, pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = pd.read_csv(sys.argv[1]); plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(figsize=(5.2, 3.8)); col = {200: '#0072B2', 1000: '#009E73', 5000: '#D55E00'}
for ct, c in col.items():
    s = d[(np.isclose(d.ctau_mm, ct)) & (d.n_ge4hits >= 10)].sort_values('Lambda_eV')
    err = np.sqrt(s.f_helical * (1 - s.f_helical) / s.n_ge4hits)
    ax.errorbar(s.Lambda_eV, s.f_helical, err, color=c, marker='o', ms=4, capsize=2, label=f'$c\\tau$ = {ct/1000:g} m' if ct >= 1000 else f'$c\\tau$ = {ct:g} mm')
if len(sys.argv) > 3:
    e = pd.read_csv(sys.argv[3])
    for ct, c in col.items():
        s = e[(np.isclose(e.ctau_mm, ct)) & (e.n_ge4hits >= 10)].sort_values('Lambda_eV')
        ax.plot(s.Lambda_eV, s.f_helical_smeared, color=c, ls='--', marker='s', ms=3, lw=1)
    ax.plot([], [], color='0.3', ls='-', label='unsmeared, RMS < 0.1 mm'); ax.plot([], [], color='0.3', ls='--', label='smeared hits, $\\chi^2$/ndf < 3')
ax.set_xscale('log'); ax.axhline(1, color='0.7', lw=.8); ax.set_xlabel(r'$\Lambda$ [eV]'); ax.set_ylabel('helical-compatible fraction ($\\geq$4 hits)')
ax.legend(frameon=False, fontsize=8); ax.set_ylim(0, 1.08); fig.tight_layout(); fig.savefig(sys.argv[2], dpi=170)
