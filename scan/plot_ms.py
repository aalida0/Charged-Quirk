"""Plots from ms_summary.csv: cross sections, terminating-track and neutral-pair yields vs mass.  usage: plot_ms.py ms_summary.csv outprefix"""
import sys, pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = pd.read_csv(sys.argv[1]); pre = sys.argv[2]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
C = {30: '#0072B2', 100: '#009E73', 300: '#E69F00', 1000: '#CC79A7'}
LUMI = 3000.  # fb^-1, HL-LHC
# (1) production cross sections
fig, ax = plt.subplots(figsize=(5.2, 3.8))
for rep, ls in (('doublet', '--'), ('triplet', '-')):
    s = d[(d.rep == rep) & (d.Lambda_eV == 100)]
    ax.plot(s.M_GeV, s.sigma_DY_fb, ls, color='#0072B2', marker='o', ms=3, label=f'{rep}: DY $Q^+Q^-$')
    ax.plot(s.M_GeV, s.sigma_Wstar_fb, ls, color='#D55E00', marker='s', ms=3, label=f'{rep}: $W^*\\to Q^\\pm Q^0$')
ax.set_yscale('log'); ax.set_xlabel('$M$ [GeV]'); ax.set_ylabel(r'$\sigma$ ($|\eta|<2.5$, LO) [fb]'); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(pre + '_xsec.png', dpi=170); plt.close(fig)
# (2) yields at HL-LHC
fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.8), sharex=True)
ax = axs[0]
MK = {30: 'o', 100: 's', 300: '^', 1000: 'D'}
for L in (30, 100, 300, 1000):
    for rep, ls in (('triplet', '-'), ('doublet', '--')):
        s_ = d[(d.rep == rep) & (d.Lambda_eV == L) & (d.sigma_termtrack_fb > 0)]
        if len(s_): ax.plot(s_.M_GeV, s_.sigma_termtrack_fb * LUMI, ls, color=C[L], marker=MK[L], ms=4, label=f'{rep}, $\\Lambda$={L} eV' if rep == 'triplet' else None, alpha=0.9)
s_ = d[(d.rep == 'doublet') & (d.Lambda_eV == 100) & (d.sigma_termtrack_fb > 0)]
ax.plot(s_.M_GeV, s_.sigma_termtrack_fb * LUMI, '--', color='0.3', lw=1, label='doublet (all $\\Lambda$ overlap; 0 events above 500 GeV)')
ax.set_title(r'terminating-track candidates ($\geq$3 tracker layers, decay in tracker)', fontsize=8)
ax.set_ylabel(r'events at 3 ab$^{-1}$ (LO, no efficiencies)'); ax.legend(frameon=False, fontsize=7)
ax = axs[1]
for rep, ls, col in (('triplet', '-', '#0072B2'), ('doublet', '--', '#D55E00')):
    g = d[d.rep == rep].groupby('M_GeV').sigma_neutralpair_MET250_fb
    ax.plot(g.mean().index, g.mean() * LUMI, ls, color=col, marker='o', ms=4, label=rep)
    ax.fill_between(g.mean().index, g.min() * LUMI, g.max() * LUMI, color=col, alpha=0.25, lw=0)
ax.set_title(r'neutral pair leaves after a decay, $p_T^{\rm pair}>250$ GeV (band: $\Lambda$ = 30...1000 eV)', fontsize=8); ax.legend(frameon=False, fontsize=8)
for ax in axs:
    ax.axhline(3., color='0.6', lw=0.8, ls=':'); ax.set_yscale('log'); ax.set_xlabel('$M$ [GeV]'); ax.set_ylim(1e-1, None)
fig.tight_layout(); fig.savefig(pre + '_yield.png', dpi=170); plt.close(fig)
