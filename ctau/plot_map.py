import numpy as np, csv, sys
sys.path.insert(0,'.')
from ctau import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
INK,SUB,GRID='#191919','#666666','#e6e6e6'
D=np.logspace(np.log10(0.12),np.log10(60),160)
cases={'doublet (Y=1/2), full strength':(1,1.0),'triplet (Y=0), full strength':(2,1.0),'doublet, coupling$^2$ x 3e-6':(1,3e-6)}
rows=[]; curves={}
for k,(kap,ex) in cases.items():
    curves[k]=np.array([ctau_cm(d,kap,ex) for d in D])
with open('ctau_map.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['dm_GeV']+[f'ctau_cm[{k}]' for k in cases])
    for i,d in enumerate(D): w.writerow([f'{d:.5g}']+[f'{curves[k][i]:.5g}' for k in cases])
cols=['#257de0','#f47230','#6a4ce2']
fig,ax=plt.subplots(figsize=(8.2,5.2),dpi=160)
ax.set_facecolor('white'); fig.patch.set_facecolor('white')
ax.axhspan(3,100,color='#f1f1f1',zorder=0)
ax.axvspan(D_LO,D_HI,color='#e9e9e9',alpha=0.6,zorder=0,lw=0)
ax.text(1.1,2e-12*1.2,'exclusive ↔ partonic\ngap (unreliable)',color=SUB,fontsize=7,ha='center',va='bottom')
ax.text(2.5,150,'inner-tracker radii (~3 cm to ~1 m)',color=SUB,fontsize=8,va='bottom')
for (k,_),c,ls in zip(cases.items(),cols,['-','-','--']):
    ax.plot(D,curves[k],color=c,lw=2,ls=ls,solid_capstyle='round',label=k)
# pure-EW benchmarks
for name,dm,kap,c in [('Higgsino-like',delta_m(.5,.5),1,cols[0]),('wino-like',delta_m(1,0),2,cols[1])]:
    y=ctau_cm(dm,kap); ax.plot([dm],[y],'o',ms=8,mfc=c,mec='white',mew=2,zorder=5)
    ax.annotate(f'{dm*1e3:.0f} MeV, {y:.2g} cm',(dm,y),xytext=(12,14) if kap==1 else (12,26),textcoords='offset points',color=INK,fontsize=8)
# literature points (hollow)
ax.plot([10],[1e-12*c_cm],'o',ms=8,mfc='white',mec=cols[2],mew=2,zorder=5)
ax.annotate('Forsyth et al. scaling\n($\\tau$=1e-12 s at 10 GeV)',(10,1e-12*c_cm),xytext=(10,-26),textcoords='offset points',ha='left',color=INK,fontsize=8)
ax.plot([40],[3e-21*c_cm],'o',ms=8,mfc='white',mec=cols[0],mew=2,zorder=5)
ax.annotate('Harnik et al. Eq. 17\n(Δm=40 GeV)',(40,3e-21*c_cm),xytext=(-16,-34),textcoords='offset points',ha='right',color=INK,fontsize=8)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.12,60); ax.set_ylim(1e-12,1e3)
ax.set_xlabel('Δm = m(Q$^+$) − m(Q$^0$)  [GeV]',color=SUB); ax.set_ylabel('c$\\tau$(Q$^+$ → Q$^0$ + X)  [cm]',color=SUB)
ax.grid(True,which='major',color=GRID,lw=0.8); ax.set_axisbelow(True)
for s in ['top','right']: ax.spines[s].set_visible(False)
for s in ['left','bottom']: ax.spines[s].set_color('#bbbbbb')
ax.tick_params(colors=SUB,labelsize=8)
leg=ax.legend(loc='lower left',frameon=False,fontsize=8,labelcolor=INK)
ax.set_title('Decay length of the charged quirk component vs mass splitting',color=INK,fontsize=11,loc='left')
fig.text(0.01,0.01,'Static-source limit, vector-like W* coupling. e/mu/tau + pi + pi-pi0 (GS rho) below 0.8 GeV; partonic above 1.5 GeV;\nlog-interpolated in between (unreliable). Finite-M kinematics changes ctau by <1% for M>60 GeV; the finite-M change of Delta m is not shown here (see ctau_vs_M.png).',color=SUB,fontsize=6.5)
fig.tight_layout(rect=(0,0.05,1,1)); fig.savefig('ctau_map.png')
print('ok', {k:[f'{curves[k][np.argmin(abs(D-x))]:.2e}' for x in (0.17,0.36,1,10)] for k in cases})
