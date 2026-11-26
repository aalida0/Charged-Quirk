"""Charged->neutral decay length of an EW-multiplet quirk, Q+ -> Q0 + X, as a function of Delta m.

Vector-like heavy current (Dirac multiplet), one-loop gauge mass splitting.  Channels: e, mu, tau (massive
leptons), pi, pi-pi0 (Gounaris-Sakurai pion form factor), partonic ud/cs above the resonance region.
Units: GeV.  Couplings are 'full strength' (kappa=1 doublet, kappa=2 triplet); extra suppression via `extra`.

Finite-M: the static-source (M >> Delta) rate is multiplied by the exact-kinematics ratio C(Delta, M), obtained
from the spin-summed heavy-fermion tensor H_{mu nu}=Tr[(p0+M0) g_mu (p+ + M+) g_nu] contracted with the
transverse (spin-1) or longitudinal (spin-0, pion pole) hadronic tensor; see docs in `rate_ratio_*`.
"""
import numpy as np
from scipy.integrate import quad

GF = 1.1663788e-5; hbar = 6.582119569e-25; c_cm = 2.99792458e10
mW, mZ = 80.377, 91.1876; cw2 = (mW / mZ) ** 2; sw2 = 1 - cw2; alpha = 1 / 127.95; a2 = alpha / sw2  # on-shell scheme
Vud2, Vcs2 = 0.9497, 0.9475
mpi, fpi = 0.13957, 0.1304      # charged pion mass; f_pi ~ 130 MeV convention
mrho, grho = 0.7755, 0.1462     # rho(770) for the GS form factor
ml = {'e': 0.000511, 'mu': 0.10566, 'tau': 1.77686}


# ----------------------------------------------------------------------------- mass splitting
def delta_m(T, Y):
    """M_Q+ - M_Q0, asymptotic M >> mW: dM_i = -(a2/2)[mW(T(T+1)-T3^2) + (mZ/cw2)(T3-Q sw2)^2]."""
    def dM(T3):
        Q = T3 + Y
        return -(a2 / 2) * (mW * (T * (T + 1) - T3 ** 2) + (mZ / cw2) * (T3 - Q * sw2) ** 2)
    return dM(0.5 if T == 0.5 else 1) - dM(-0.5 if T == 0.5 else 0)


def _fvec(r):
    """f(r) = -(1/2pi) int_0^1 (2-x) ln((1-x)^2 + x r^2) dx : massive-vector self-energy on shell, mass M=1."""
    return -quad(lambda x: (2 - x) * np.log((1 - x) ** 2 + x * r * r), 0, 1, limit=200, points=[1 - 1e-6])[0] / (2 * np.pi)


def delta_m_M(T, Y, M):
    """One-loop splitting at finite M (Feynman gauge, Dirac multiplet, internal fermion mass = M): the Lambda- and
    scale-dependent parts are member-independent (sum of couplings = Casimir) and cancel in the difference.
    dM_i = M * sum_V alpha_V C_V(i) f(m_V/M).  Reduces to delta_m(T,Y) for M >> mW."""
    def dM(T3):
        Q = T3 + Y
        cW = a2 * (T * (T + 1) - T3 ** 2); cZ = a2 * (T3 - Q * sw2) ** 2 / cw2; cg = a2 * sw2 * Q ** 2
        return M * (cW * _fvec(mW / M) + cZ * _fvec(mZ / M) + cg * _fvec(0.))
    hi, lo = (0.5, -0.5) if T == 0.5 else (1, 0)
    return dM(hi) - dM(lo)


# ----------------------------------------------------------------------------- pion form factor and R_pipi
def _k(s):  # pion momentum in the dipion rest frame
    return 0.5 * np.sqrt(max(s - 4 * mpi ** 2, 0.))


def _h(s):
    k = _k(s); rs = np.sqrt(s)
    return (2 / np.pi) * (k / rs) * np.log((rs + 2 * k) / (2 * mpi)) if k > 0 else 0.


_k0 = _k(mrho ** 2); _h0 = _h(mrho ** 2)
_hp0 = _h0 * (1 / (8 * _k0 ** 2) - 1 / (2 * mrho ** 2)) + 1 / (2 * np.pi * mrho ** 2)
_d = (3 / np.pi) * (mpi ** 2 / _k0 ** 2) * np.log((mrho + 2 * _k0) / (2 * mpi)) + mrho / (2 * np.pi * _k0) \
    - mpi ** 2 * mrho / (np.pi * _k0 ** 3)


def FF2(s):
    """|F_pi(s)|^2, Gounaris-Sakurai rho(770) (pure rho; no omega mixing, no rho').  Defined for s >= 4 m_pi^2 only
    (R_pipi vanishes below; below threshold the value at threshold is returned and is never used in a rate)."""
    s = max(s, 4 * mpi ** 2 + 1e-12)
    k = _k(s); rs = np.sqrt(s)
    gam = grho * (mrho / rs) * (k / _k0) ** 3
    f = grho * mrho ** 2 / _k0 ** 3 * (k ** 2 * (_h(s) - _h0) + (mrho ** 2 - s) * _k0 ** 2 * _hp0)
    den = (mrho ** 2 - s + f) ** 2 + (mrho * gam) ** 2
    return (mrho ** 2 + _d * mrho * grho) ** 2 / den


def Rpipi(s):
    """Spectral ratio to one massless V-A lepton channel:  R = 3 V_ud^2 v_1(s),  v_1 = beta^3 |F|^2 / 12."""
    if s <= 4 * mpi ** 2:
        return 0.
    b = np.sqrt(1 - 4 * mpi ** 2 / s)
    return Vud2 * b ** 3 * FF2(s) / 4.


def tau_pipi_check():
    """B(tau -> pi pi0 nu)/B(tau -> e nu nu) = int dx 2(1-x)^2(1+2x) R_pipi(x m_tau^2)  (expect ~1.43 incl. rho', ~1.3-1.4 pure rho)."""
    mt2 = ml['tau'] ** 2
    return quad(lambda x: 2 * (1 - x) ** 2 * (1 + 2 * x) * Rpipi(x * mt2), 4 * mpi ** 2 / mt2, 1, limit=200)[0]


# ----------------------------------------------------------------------------- finite-M kinematics
def _gT(s, D, M0):
    """|p*| T / M+^2 for a spin-1 (transverse) hadronic tensor, exact kinematics (stable form).  Static limit: 8 (D^2-s)^{3/2}."""
    Mp = M0 + D; Dd = D * (Mp + M0)          # M+^2 - M0^2
    lam = (Dd - s) ** 2 - 4 * M0 ** 2 * s
    if lam <= 0: return 0.
    T = 4 * ((Dd ** 2 - s ** 2) / 2 + s * ((D ** 2 - s) / 2 - 2 * M0 * Mp))
    return np.sqrt(lam) / (2 * Mp) * T / Mp ** 2


def _gL(s, D, M0):
    """same for the spin-0 (q^mu q^nu, pion pole) tensor.  Static limit: 8 D^2 sqrt(D^2-s)."""
    Mp = M0 + D; Dd = D * (Mp + M0)
    lam = (Dd - s) ** 2 - 4 * M0 ** 2 * s
    if lam <= 0: return 0.
    T = 4 * ((Dd ** 2 - s ** 2) / 2 - s * (D ** 2 - s) / 2)
    return np.sqrt(lam) / (2 * Mp) * T / Mp ** 2


def C_T(D, M0, R=lambda s: 1.0, smin=0.):
    """exact/static ratio for a continuum channel with spectral ratio R(s) over s in [smin, D^2]."""
    if M0 is None or D <= np.sqrt(smin): return 1.0
    num = quad(lambda s: R(s) * _gT(s, D, M0), smin, D * D, limit=200)[0]
    den = quad(lambda s: R(s) * 8 * (D * D - s) ** 1.5, smin, D * D, limit=200)[0]
    return num / den if den > 0 else 1.0


def C_L(D, M0, s=mpi ** 2):
    if M0 is None or D <= np.sqrt(s): return 1.0
    return _gL(s, D, M0) / (8 * D * D * np.sqrt(D * D - s))


# ----------------------------------------------------------------------------- widths
def F2(D, ma, mb):
    """30/D^5 * int E_a p_a E_b p_b dE_a, E_b=D-E_a ; ->1 for massless"""
    if D <= ma + mb: return 0.
    f = lambda E: (E * np.sqrt(max(E * E - ma * ma, 0))) * ((D - E) * np.sqrt(max((D - E) ** 2 - mb * mb, 0)))
    return 30 / D ** 5 * quad(f, ma, D - mb, limit=200)[0]


def g_pipi(D, kappa=1.0, M0=None):
    """Gamma(Q+ -> Q0 pi+ pi0): (5/2) int ds (D^2-s)^{3/2} R_pipi(s) x Gamma_lept(massless)/D^5  [x C_T]."""
    if D <= 2 * mpi: return 0.
    base = kappa * GF ** 2 / (15 * np.pi ** 3)
    I = quad(lambda s: (D * D - s) ** 1.5 * Rpipi(s), 4 * mpi ** 2, D * D, limit=200)[0]
    return base * 2.5 * I * C_T(D, M0, Rpipi, 4 * mpi ** 2)


def channel_widths(D, kappa=1.0, M0=None):
    """exclusive partial widths [GeV] by channel: {'e','mu','tau','pi','pipi'}; valid for D <~ 0.8 GeV."""
    base = kappa * GF ** 2 * D ** 5 / (15 * np.pi ** 3); cT = C_T(D, M0)
    out = {l: base * F2(D, ml[l], 0.) * cT for l in ml}
    out['pi'] = kappa * GF ** 2 * Vud2 * fpi ** 2 * D ** 3 * np.sqrt(1 - (mpi / D) ** 2) / np.pi * C_L(D, M0) if D > mpi else 0.
    out['pipi'] = g_pipi(D, kappa, M0)
    return out


def widths(D, kappa=1.0, M0=None, pipi=True):
    """(exclusive, partonic) total widths in GeV.  exclusive = l + pi (+ pi pi0).  M0=None: static limit."""
    base = kappa * GF ** 2 * D ** 5 / (15 * np.pi ** 3)
    cT = C_T(D, M0)
    lept = sum(base * F2(D, ml[l], 0.) for l in ml) * cT
    gpi = 0.
    if D > mpi:
        gpi = kappa * GF ** 2 * Vud2 * fpi ** 2 * D ** 3 * np.sqrt(1 - (mpi / D) ** 2) / np.pi * C_L(D, M0)
    gpp = g_pipi(D, kappa, M0) if pipi else 0.
    excl = lept + gpi + gpp
    part = lept + 3 * base * (Vud2 * F2(D, 0., 0.) + Vcs2 * F2(D, 1.27, 0.095)) * cT
    return excl, part


# exclusive channels (pi, pi pi, l) are used up to D_LO; partonic above D_HI; in between both prescriptions bracket the answer
D_LO, D_HI = 0.8, 1.5


def ctau_cm(D, kappa=1.0, extra=1.0, M0=None, pipi=True, band=False):
    """hybrid: exclusive for D<=D_LO; partonic for D>=D_HI; log-linear interpolation in between (UNRELIABLE region).
    band=True returns (central, lo, hi) with lo/hi the smaller/larger of the two prescriptions inside the gap."""
    def ct(G): return hbar / (G * extra) * c_cm
    e, p = widths(D, kappa, M0, pipi)
    if D <= D_LO: out = (ct(e),) * 3
    elif D >= D_HI: out = (ct(p),) * 3
    else:
        a = ct(widths(D_LO, kappa, M0, pipi)[0]); b = ct(widths(D_HI, kappa, M0, pipi)[1])
        t = (np.log(D) - np.log(D_LO)) / (np.log(D_HI) - np.log(D_LO))
        c = float(np.exp((1 - t) * np.log(a) + t * np.log(b)))
        out = (c, min(ct(e), ct(p)), max(ct(e), ct(p)))
    return out if band else out[0]


if __name__ == "__main__":
    print(f"tau->pipi nu ratio to B_e (GS rho only, normalisation check): {tau_pipi_check():.3f}")
    print(f"|F(m_rho^2)|^2={FF2(mrho**2):.2f}  |F(0.3)|^2={FF2(0.3**2):.3f}  |F(0.1)|^2 @thr={FF2(4*mpi**2):.3f}")
    for name, T, Y, k in [("doublet", .5, .5, 1), ("triplet", 1, 0., 2)]:
        d0 = delta_m(T, Y)
        print(f"{name}: asymptotic dm = {d0*1e3:.1f} MeV;  finite-M dm: " +
              "  ".join(f"M={M:g}: {delta_m_M(T,Y,M)*1e3:.1f}" for M in (50, 100, 200, 500, 1000, 1e4)))
        e0, _ = widths(d0, k)
        e1, _ = widths(d0, k, M0=500.)
        e2, _ = widths(d0, k, M0=100.)
        print(f"   ctau static={hbar/e0*c_cm:.3f} cm ; M0=500: {hbar/e1*c_cm:.3f} cm ; M0=100: {hbar/e2*c_cm:.3f} cm")
