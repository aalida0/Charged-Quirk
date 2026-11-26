"""Demo plots from a quirksim triplet decay run (see quirkdecay/macros/triplet.mac).
usage: python3 plot_triplet_demo.py trip.root out_prefix M_GeV ctau_mm"""
import sys, numpy as np, uproot, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
f, out, M, ctau = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
a = uproot.open(f)['events'].arrays(['dec_x','dec_y','dec_z','dec_ppx','dec_ppy','dec_ppz','dec_pe',
        'dec_dpx','dec_dpy','dec_dpz','dec_de','dp_pdg','dp_px','dp_py','dp_pz','dp_e'], library='np')
cat = lambda k: np.concatenate([np.asarray(x, float) for x in a[k]])
x, y, z = cat('dec_x'), cat('dec_y'), cat('dec_z')          # mm
ppx, ppy, ppz = cat('dec_ppx'), cat('dec_ppy'), cat('dec_ppz')  # parent momentum (GeV)
r = np.hypot(x, y); L = np.sqrt(x**2+y**2+z**2)
pp = np.sqrt(ppx**2+ppy**2+ppz**2)/1e3   # ntuple momenta are in MeV; pp in GeV
pdg = cat('dp_pdg'); px, py, pz = cat('dp_px'), cat('dp_py'), cat('dp_pz')
pt = np.hypot(px, py); pmag = np.sqrt(px**2+py**2+pz**2)
ispi = np.abs(pdg) == 211
# kink: angle between parent and neutral-daughter momenta (dec_d* = neutral quirk, MeV)
dpx, dpy, dpz = cat('dec_dpx'), cat('dec_dpy'), cat('dec_dpz')
print(f'decays {len(L)}, mean 3D decay distance {L.mean():.1f} mm, frac r<1.2 m {np.mean(r<1200):.3f}')
fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
c = '#2a6f97'
bins = np.linspace(0, 1000, 41)
ax[0].hist(L, bins, color=c, alpha=.85, density=True)
lam = np.mean(pp/M*ctau)
xs = np.linspace(0, 1000, 200); ax[0].plot(xs, np.exp(-xs/lam)/lam, 'k--', lw=1)
ax[0].set(xlabel='decay distance from vertex [mm]', ylabel='density [1/mm]', title='Q$^+$ decay point')
ax[0].text(.55, .8, f'$\\langle\\beta\\gamma\\rangle c\\tau$ = {lam:.0f} mm\n(dashed: straight-line\nexponential)', transform=ax[0].transAxes, fontsize=8)
ax[1].hist(pt[ispi], np.linspace(0, 400, 41), color=c, alpha=.85)
ax[1].set(xlabel='$\\pi^\\pm$ lab $p_T$ [MeV]', ylabel='decays', title='soft decay pion')
ang = np.arccos(np.clip((ppx*dpx+ppy*dpy+ppz*dpz)/(np.sqrt(ppx**2+ppy**2+ppz**2)*np.sqrt(dpx**2+dpy**2+dpz**2)), -1, 1))*1e3
ax[2].hist(ang, 40, color=c, alpha=.85)
ax[2].set(xlabel='kink angle Q$^+\\to$Q$^0$ [mrad]', ylabel='decays', title='track kink')
for s in ('top','right'):
    for q in ax: q.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(out+'.png', dpi=200)
print(f'pion fraction {ispi.sum()/len(L):.3f}; median pi pT {np.median(pt[ispi]):.0f} MeV; kink median {np.median(ang):.3f} mrad, 95% {np.percentile(ang,95):.3f}')
