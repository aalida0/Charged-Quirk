"""LO re-weighting of the Pythia gamma*/Z -> 17 17bar events (SM 4th-generation lepton couplings, V-A) to a vector-like multiplet member.
usage: reweight_dy.py <dyf.hepmc3> <out.npz> <M_GeV>
For every event: weight = sigma_hat(q qbar -> Q+ Q-; ŝ, target couplings) / sigma_hat(q qbar -> L- L+; ŝ, Pythia couplings), with
  sigma_hat ∝ beta * sum_{lambda=L,R} [ (3-beta^2)/2 |C_V,lambda|^2 + beta^2 |C_A,lambda|^2 ],
  C_V,lambda = Q_q Q_F / ŝ + g^q_lambda g^F_V /(sw2 cw2 (ŝ - mZ^2 + i mZ GZ)),  C_A,lambda = g^q_lambda g^F_A /(sw2 cw2 (ŝ - mZ^2 + i mZ GZ))
(g = T3 - Q sw2 couplings; the Z coupling is e/(sw cw) g, hence the 1/(sw2 cw2); massless quarks).  The ŝ, the incoming flavours and the event numbers are read from GenPdfInfo.
Targets: 'doublet' (T=1/2, Y=1/2, Q+: g_V = 1/2 - sw2), 'triplet' (T=1, Y=0, Q+: g_V = 1 - sw2), vector-like (g_A = 0);
reference: charged lepton (g_L = -1/2 + sw2, g_R = sw2, i.e. g_V = (g_L+g_R)/2, g_A = (g_R-g_L)/2 up to sign)."""
import sys, re, numpy as np
fin, fout, M = sys.argv[1], sys.argv[2], float(sys.argv[3])
sw2, mZ, GZ = 0.2312, 91.1876, 2.4952
SQRTS = 13000.
q = {1: (-1/3, -.5), 2: (2/3, .5), 3: (-1/3, -.5), 4: (2/3, .5), 5: (-1/3, -.5)}   # (charge, T3L)
def sig(flav, shat, QF, gV, gA):
    Qq, T3 = q[flav]; gqL, gqR = T3 - Qq * sw2, -Qq * sw2
    b = np.sqrt(max(1 - 4 * M * M / shat, 1e-12)); PZ = 1 / (sw2 * (1 - sw2) * (shat - mZ ** 2 + 1j * mZ * GZ))
    tot = 0.
    for gq in (gqL, gqR):
        CV = Qq * QF / shat + gq * gV * PZ; CA = gq * gA * PZ
        tot += (3 - b * b) / 2 * abs(CV) ** 2 + b * b * abs(CA) ** 2
    return b * tot
ref = dict(QF=-1., gV=((-.5 + sw2) + sw2) / 2, gA=(sw2 - (-.5 + sw2)) / 2)
tgt = {'doublet': dict(QF=-1., gV=-(0.5 - sw2), gA=0.), 'triplet': dict(QF=-1., gV=-(1 - sw2), gA=0.)}   # antiparticle conventions: overall sign flips both QF and gV
ev, w = [], {k: [] for k in tgt}
cur = None
def flush():
    if cur is None or cur.get('id1') is None: return
    shat = cur['x1'] * cur['x2'] * SQRTS ** 2; fl = abs(cur['id1'])
    if fl not in q: return
    ev.append(cur['n'])
    for k, t in tgt.items(): w[k].append(sig(fl, shat, **t) / sig(fl, shat, **ref))
for l in open(fin):
    if l.startswith('E '):
        flush(); cur = {'n': int(l.split()[1])}
    elif l.startswith('A 0 GenPdfInfo'):
        f = l.split(); cur['id1'], cur['id2'], cur['x1'], cur['x2'] = int(f[3]), int(f[4]), float(f[5]), float(f[6])
flush()
np.savez(fout, evnum=np.array(ev), **{k: np.array(v) for k, v in w.items()})
for k in tgt: print(f'{k}: <w> = {np.mean(w[k]):.3f}  (n = {len(ev)})')
