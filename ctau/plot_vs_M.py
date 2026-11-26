"""c*tau(Q+) and Delta m for the pure-EW doublet (Y=1/2) and triplet (Y=0) quirk as a function of the neutral mass M:
one-loop Delta m(M) (Ibe-Matsumoto-Sato form, checked) and exact-kinematics widths."""
import numpy as np, csv
from ctau import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M=np.logspace(np.log10(50),np.log10(5000),40)
cols={'doublet':'#257de0','triplet':'#f47230'}; rep={'doublet':(.5,.5,1),'triplet':(1,0.,2)}
out={}
for n,(T,Y,k) in rep.items():
    dm=np.array([delta_m_M(T,Y,m) for m in M]); ct=np.array([hbar/widths(d,k,M0=m)[0]*c_cm for d,m in zip(dm,M)])
    out[n]=(dm,ct)
with open('ctau_vs_M.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['M_GeV']+[f'{n}_{q}' for n in rep for q in ('dm_MeV','ctau_cm')])
    for i,m in enumerate(M): w.writerow([f'{m:.5g}']+[f'{v:.5g}' for n in rep for v in (out[n][0][i]*1e3,out[n][1][i])])
INK,SUB,GRID='#191919','#666666','#e6e6e6'
fig,ax=plt.subplots(figsize=(6.6,4.4),dpi=160)
for n,c in cols.items():
    ax.plot(M,out[n][1],color=c,lw=2,label=f'{n}  (Δm: {out[n][0][0]*1e3:.0f} → {out[n][0][-1]*1e3:.0f} MeV)')
    ax.axhline(hbar/widths(delta_m(*rep[n][:2]),rep[n][2])[0]*c_cm,color=c,lw=1,ls=':')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('M(Q$^0$) [GeV]',color=SUB); ax.set_ylabel('c$\\tau$(Q$^+$) [cm]',color=SUB)
ax.grid(True,which='major',color=GRID,lw=.8); ax.set_axisbelow(True)
for s_ in ['top','right']: ax.spines[s_].set_visible(False)
ax.tick_params(colors=SUB,labelsize=8); ax.legend(frameon=False,fontsize=8,labelcolor=INK)
ax.set_title('Pure-EW multiplet: decay length vs mass (dotted: M → ∞ value)',loc='left',fontsize=10,color=INK)
fig.tight_layout(); fig.savefig('ctau_vs_M.png'); print('ok')
