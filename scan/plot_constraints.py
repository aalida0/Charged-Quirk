"""Combined (Lambda, M) constraint map for N_ic = 2, 3 (referee v3 P0-E). Inputs: atlas_dt_limits_intended.csv (map-based recast, 'string' rows, raw..calibrated band),
atlas2026_scaling.txt values (scaling ESTIMATE for helical tracks only, Lambda<=100 eV), pilot_joint2*.csv (non-helical candidate yields, chi2/ndf>9 flag).
usage: plot_constraints.py <data/massscan> <out.png>"""
import sys, glob, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
dd, out = sys.argv[1:3]
lim = pd.read_csv(f'{dd}/atlas_dt_limits_intended.csv'); lim = lim[lim.treatment == 'string']
pj = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f'{dd}/pilot_joint2*.csv'))]).drop_duplicates(['M_GeV', 'Lambda_eV'], keep='last')
est = {2: (872, 986), 3: (936, 1054)}                      # atlas2026_scaling.py: obs-scaled, exp-scaled (helical only)
col = dict(dt='#0072B2', est='#D55E00', pj='#009E73')
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
for ax, nic in zip(axs, (2, 3)):
    l = lim[lim.N_IC == nic].sort_values('Lambda_eV'); L = l.Lambda_eV.values
    ax.fill_between(L, 250, l.limit_raw_GeV, color=col['dt'], alpha=.12, lw=0); ax.fill_between(L, l.limit_raw_GeV, l.limit_calibrated_GeV, color=col['dt'], alpha=.35, lw=0)
    ax.plot(L, l.limit_raw_GeV, color=col['dt'], lw=1.2, ls='--'); ax.plot(L, l.limit_calibrated_GeV, color=col['dt'], lw=1.2)
    ax.annotate('', xy=(1000, 300), xytext=(1000, 380), arrowprops=dict(arrowstyle='-|>', color=col['dt'], lw=1.2))
    ax.fill_between([1, 100], *est[nic], color=col['est'], alpha=.25, lw=0, hatch='///', edgecolor=col['est'])
    for _, r in pj.iterrows():
        c = pj[(pj.M_GeV == r.M_GeV) & (pj.Lambda_eV <= 100)]                      # pooled helical controls (Lambda = 1, 30, 100 eV)
        ctrl = c['nEndpoints_chi2gt9'].sum() / max(c['nEndpoints_L4'].sum(), 1)    # chi2>9 false-positive fraction of the control
        n9 = (r['sigma_C9_fb'] - ctrl * r['sigma_L4_fb']) * 3000 * nic      # excess events at 3/ab, before any reco efficiency, times N_ic
        if n9 > 5 and r.Lambda_eV >= 300: ax.scatter(r.Lambda_eV, r.M_GeV, s=7 * np.log10(n9) ** 2, facecolors='none', edgecolors=col['pj'], lw=1.1, zorder=5)
    ax.axhline(300, color='gray', lw=.6, ls=':')
    ax.set_xscale('log'); ax.set_xlim(1, 2000); ax.set_ylim(230, 1150); ax.set_xlabel(r'$\Lambda$ (eV)'); ax.set_title(rf'$N_{{\rm ic}}={nic}$', fontsize=9)
    ax.tick_params(labelsize=8)
axs[0].set_ylabel(r'$M$ (GeV), Dirac triplet')
h = [plt.Line2D([], [], color=col['dt'], lw=1.2), plt.Rectangle((0, 0), 1, 1, fc=col['dt'], alpha=.35), plt.Rectangle((0, 0), 1, 1, fc=col['est'], alpha=.25, hatch='///', ec=col['est']),
     plt.Line2D([], [], marker='o', ls='', mfc='none', mec=col['pj'])]
fig.legend(h, ['2018-19 recast, calibrated (dashed: raw)', 'excluded (raw to calibrated)', '2026 search, scaling estimate, helical only', 'excess non-helical endpoints over the string-free control (area grows with log of the 3 ab$^{-1}$ count)'],
           loc='lower center', ncol=2, fontsize=6.5, frameon=False)
fig.tight_layout(rect=(0, .14, 1, 1)); fig.savefig(out, dpi=200)
