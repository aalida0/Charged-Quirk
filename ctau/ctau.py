"""Charged->neutral decay length of an EW-multiplet quirk, Q+ -> Q0 + X, as a function of Delta m.
Static-source (M >> Delta) limit, vector-like heavy current, one-loop gauge mass splitting.
Units: GeV. All couplings are 'full strength' (kappa=1 doublet, kappa=2 triplet); extra suppression via kappa_extra.
"""
import numpy as np
from scipy.integrate import quad

GF=1.1663788e-5; hbar=6.582119569e-25; c_cm=2.99792458e10
mW,mZ=80.377,91.1876; cw2=(mW/mZ)**2; sw2=1-cw2; alpha=1/127.95; a2=alpha/sw2  # on-shell scheme
Vud2,Vcs2=0.9497,0.9475
mpi,fpi=0.13957,0.1304      # f_pi ~130 MeV convention
ml={'e':0.000511,'mu':0.10566,'tau':1.77686}

def delta_m(T,Y):
    """M_Q - M_0 (asymptotic M>>mW), from massive W/Z Yukawa self-energy: dM_i = -(a2/2)[mW(T(T+1)-T3^2)+(mZ/cw2)(T3-Q sw2)^2]"""
    def dM(T3):
        Q=T3+Y
        return -(a2/2)*(mW*(T*(T+1)-T3**2)+(mZ/cw2)*(T3-Q*sw2)**2)
    return dM(T3=Y*0+ (0.5 if T==0.5 else 1)) - dM(T3=(-0.5 if T==0.5 else 0))   # charged minus neutral

def F2(D,ma,mb):
    """30/D^5 * int E_a p_a E_b p_b dE_a, E_b=D-E_a ; ->1 for massless"""
    if D<=ma+mb: return 0.
    f=lambda E:(E*np.sqrt(max(E*E-ma*ma,0)))*((D-E)*np.sqrt(max((D-E)**2-mb*mb,0)))
    return 30/D**5*quad(f,ma,D-mb,limit=200)[0]

def widths(D,kappa=1.0):
    """returns (exclusive-ish, partonic) total widths in GeV for coupling factor kappa (1=doublet,2=triplet)"""
    base=kappa*GF**2*D**5/(15*np.pi**3)
    lept=sum(base*F2(D,ml[l],0.) for l in ml)
    # exclusive hadronic: pi only (valid D<~0.5)
    gpi=0.
    if D>mpi:
        gpi=kappa*GF**2*Vud2*fpi**2*D**3*np.sqrt(1-(mpi/D)**2)/np.pi
    excl=lept+gpi
    # partonic hadronic: 3*Vud^2 (m=0) + 3*Vcs^2 (m_c=1.27,m_s=0.095)
    part=lept+3*base*(Vud2*F2(D,0.,0.)+Vcs2*F2(D,1.27,0.095))
    return excl,part

def ctau_cm(D,kappa=1.0,extra=1.0):
    """hybrid: exclusive (pi,e,mu) for D<=0.5; partonic for D>=1.5; log-linear interpolation of ctau in between (UNRELIABLE region)"""
    def ct(G): return hbar/(G*extra)*c_cm
    if D<=0.5: return ct(widths(D,kappa)[0])
    if D>=1.5: return ct(widths(D,kappa)[1])
    a=ct(widths(0.5,kappa)[0]); b=ct(widths(1.5,kappa)[1])
    t=(np.log(D)-np.log(0.5))/(np.log(1.5)-np.log(0.5))
    return float(np.exp((1-t)*np.log(a)+t*np.log(b)))

if __name__=="__main__":
    dd=delta_m(0.5,0.5); dt=delta_m(1,0.)
    print(f"Delta m doublet Y=1/2: {dd*1e3:.0f} MeV ; triplet Y=0: {dt*1e3:.0f} MeV  (asymptotic, M>>mW)")
    for name,D,k in [("doublet",dd,1),("triplet",dt,2)]:
        e,p=widths(D,k)
        print(f"{name}: width(excl)={e:.3e} GeV  tau={hbar/e:.3e} s  ctau={hbar/e*c_cm:.3f} cm")
    # cross-checks vs Harnik Eq.17 and Forsyth scaling at D=10 GeV
    e,p=widths(10.,1.)
    print(f"D=10 GeV full-strength doublet (partonic): tau={hbar/p:.2e} s ; Harnik Eq17 (s_theta=0.717): {3e-21*(40/10)**5:.2e} s ; Forsyth: 1e-12 s")
    for D in [0.2,0.35,0.5,1,2,5,10,40]:
        print(f"D={D:5.2f}  ctau(kappa=1)={ctau_cm(D,1):.3e} cm   ctau(kappa=2)={ctau_cm(D,2):.3e} cm")
