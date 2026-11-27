"""Generator-level MET and anti-kT(0.4) jets per HepMC3 event (Pythia hadron level), for a monojet-style selection of the neutral pair.
usage: jetmet_gen.py <in.hepmc3> <out.npz>
Visible = status-1 particles that are not neutrinos and not 17/18 (quirks), |eta| < 4.7 for MET; jets from visible particles with pT > 1 GeV, |eta| < 4.7.
Outputs per event: evnum, met (GeV), jet pT/eta/phi (up to 4 leading jets with pT > 30 GeV, eta < 4.7), n_jets30."""
import sys, math, numpy as np
INV = {12, 14, 16, 17, 18}
def antikt(pt, eta, phi, R=0.4, ptmin=30.):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    P = np.c_[px, py, pz, np.sqrt(px ** 2 + py ** 2 + pz ** 2)]
    def kin(P):
        ptv = np.hypot(P[:, 0], P[:, 1]); y = 0.5 * np.log((P[:, 3] + P[:, 2]) / np.maximum(P[:, 3] - P[:, 2], 1e-9)); ph = np.arctan2(P[:, 1], P[:, 0]); return ptv, y, ph
    jets = []
    while len(P):
        ptv, y, ph = kin(P); n = len(P)
        dphi = np.abs(ph[:, None] - ph[None, :]); dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
        dR2 = (y[:, None] - y[None, :]) ** 2 + dphi ** 2
        ik = 1. / np.maximum(ptv, 1e-9) ** 2
        dij = np.minimum(ik[:, None], ik[None, :]) * dR2 / R ** 2; np.fill_diagonal(dij, np.inf)
        iB = ik
        a = np.unravel_index(np.argmin(dij), dij.shape); dmin = dij[a]
        if dmin < iB.min():
            i, j = a; P[i] = P[i] + P[j]; P = np.delete(P, j, axis=0)
        else:
            k = int(np.argmin(iB)); jets.append(P[k].copy()); P = np.delete(P, k, axis=0)
    out = []
    for J in jets:
        ptj = math.hypot(J[0], J[1])
        if ptj > ptmin: out.append((ptj, math.asinh(J[2] / ptj), math.atan2(J[1], J[0])))
    return sorted(out, reverse=True)
ev, met, jets, nj = [], [], [], []
cur = None; parts = []
def flush():
    if cur is None: return
    a = np.array(parts, float).reshape(-1, 4)  # px py pz e
    if len(a):
        pt = np.hypot(a[:, 0], a[:, 1]); eta = np.arcsinh(a[:, 2] / np.maximum(pt, 1e-9))
        m = (np.abs(eta) < 4.7)
        met.append(math.hypot(a[m, 0].sum(), a[m, 1].sum()))
        s = m & (pt > 1.)
        js = antikt(pt[s], eta[s], np.arctan2(a[s, 1], a[s, 0])) if s.sum() else []
    else: met.append(0.); js = []
    ev.append(cur); nj.append(len(js)); jets.append((js + [(0., 0., 0.)] * 4)[:4])
for l in open(sys.argv[1]):
    if l.startswith('E '):
        flush(); cur = int(l.split()[1]); parts = []
    elif l.startswith('P ') and l.split()[-1] == '1':
        f = l.split(); pdg = abs(int(f[3]))
        if pdg in INV: continue
        parts.append(tuple(map(float, f[4:8])))
flush()
np.savez(sys.argv[2], evnum=np.array(ev), met=np.array(met), jets=np.array(jets), nj30=np.array(nj))
print('events', len(ev), 'median MET', np.median(met))
