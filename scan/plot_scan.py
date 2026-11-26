"""Heat maps of the Lambda x c*tau scan + analytic string-loss contours. usage: plot_scan.py scan.csv out_prefix"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
df = pd.read_csv(sys.argv[1]); out = sys.argv[2]
Ls = np.sort(df.Lambda_eV.unique()); Cs = np.sort(df.ctau_mm.unique())
def grid(col):
    g = df.pivot(index='Lambda_eV', columns='ctau_mm', values=col).reindex(index=Ls, columns=Cs)
    return g.values
HBARC = 1.973e-13  # GeV mm
def tension(L_eV): return (L_eV*1e-9)**2 / HBARC   # GeV/mm
panels = [('f_tracker', 'decays inside the tracker', 'Blues', (0, 1)),
          ('f_nl_ge3', r'decays after $\geq 3$ tracker layers hit', 'Blues', (0, 1)),
          ('dke_median', r'median $\Delta E_{kin}/E_{kin}$ lost before the decay', 'Oranges', (0, 1)),
          ('f_survive', 'survive to the calorimeters (HSCP-like)', 'Greens', (0, 1))]
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
X, Y = np.arange(len(Cs)), np.arange(len(Ls))
for ax, (col, title, cm, (lo, hi)) in zip(axs.ravel(), panels):
    g = grid(col); im = ax.imshow(g, origin='lower', cmap=cm, vmin=lo, vmax=hi, aspect='auto')
    for i in range(len(Ls)):
        for j in range(len(Cs)):
            v = g[i, j]
            if np.isfinite(v): ax.text(j, i, '0' if v == 0 else (f'{v:.2f}' if v >= 0.01 else f'{v:.0e}'.replace('e-0', 'e-')), ha='center', va='center',
                                      fontsize=7, color='white' if v > 0.6 else 'black')
    ax.set_xticks(X); ax.set_xticklabels([f'{c:g}' for c in Cs]); ax.set_yticks(Y); ax.set_yticklabels([f'{l:g}' for l in Ls])
    ax.set_xlabel(r'$c\tau$ [mm]'); ax.set_ylabel(r'$\Lambda$ [eV]'); ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=.046)
plt.tight_layout(); plt.savefig(out + '_maps.png', dpi=180)

# string-loss check: simulated median energy loss vs T * decay distance
fig, ax = plt.subplots(figsize=(5.2, 4.2))
d = df[(df.f_decayed > 0.5) & (df.Ld_med_mm > 0)].copy()
d['ana'] = tension(d.Lambda_eV) * d.Ld_med_mm
sc = ax.scatter(d.ana, d.dE_med_GeV, c=np.log10(d.Lambda_eV), cmap='viridis', s=28)
lim = np.array([1e-6, 1e3]); ax.plot(lim, lim, 'k--', lw=1, label='1'); ax.plot(lim, lim*0.72, 'k:', lw=1, label='0.72')
ax.legend(title=r'$\Delta E / (T L_d)$', fontsize=8, loc='upper left')
ax.set(xscale='log', yscale='log', xlim=(1e-5, 1e3), ylim=(1e-5, 1e3), xlabel=r'$T\,L_d$, $T=\Lambda^2/\hbar c$, median decay distance $L_d$ [GeV]',
       ylabel='median simulated energy loss [GeV]'); plt.colorbar(sc, label=r'$\log_{10}\Lambda/$eV')
plt.tight_layout(); plt.savefig(out + '_loss_check.png', dpi=180)
ok = d[d.dE_med_GeV / d.ke0_med_GeV < 0.1]
print('ratio dE/(T Ld) for dE/KE<0.1: median %.3f, min %.3f, max %.3f (n=%d)' % ((ok.dE_med_GeV/ok.ana).median(), (ok.dE_med_GeV/ok.ana).min(), (ok.dE_med_GeV/ok.ana).max(), len(ok)))
print(d[['Lambda_eV', 'ctau_mm', 'Ld_med_mm', 'ke0_med_GeV', 'dE_med_GeV', 'ana']].assign(ratio=d.dE_med_GeV / d.ana).round(4).to_string())
