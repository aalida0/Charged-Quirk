"""LO re-weighting of the Pythia gamma*/Z -> 17 17bar events to the NEUTRAL member Q0 of the vector-like doublet (Z only, g_V = T3 = 1/2, Q = 0).
Same method as reweight_dy.py (copy with an additional target); usage: reweight_dy_n0.py <rec_dy.hepmc3> <out.npz> <M_GeV>"""
import sys, numpy as np
fin, fout, M = sys.argv[1], sys.argv[2], float(sys.argv[3])
sw2, mZ, GZ = 0.2312, 91.1876, 2.4952; SQRTS = 13000.
q = {1: (-1/3, -.5), 2: (2/3, .5), 3: (-1/3, -.5), 4: (2/3, .5), 5: (-1/3, -.5)}
def sig(flav, shat, QF, gV, gA):
    Qq, T3 = q[flav]; gqL, gqR = T3 - Qq * sw2, -Qq * sw2
    b = np.sqrt(max(1 - 4 * M * M / shat, 1e-12)); PZ = 1 / (sw2 * (1 - sw2) * (shat - mZ ** 2 + 1j * mZ * GZ)); tot = 0.
    for gq in (gqL, gqR):
        CV = Qq * QF / shat + gq * gV * PZ; CA = gq * gA * PZ
        tot += (3 - b * b) / 2 * abs(CV) ** 2 + b * b * abs(CA) ** 2
    return b * tot
ref = dict(QF=-1., gV=((-.5 + sw2) + sw2) / 2, gA=(sw2 - (-.5 + sw2)) / 2)
tgt = dict(QF=0., gV=0.5, gA=0.)
ev, w, cur = [], [], None
def flush():
    if cur is None or cur.get('id1') is None: return
    shat = cur['x1'] * cur['x2'] * SQRTS ** 2; fl = abs(cur['id1'])
    if fl not in q: return
    ev.append(cur['n']); w.append(sig(fl, shat, **tgt) / sig(fl, shat, **ref))
for l in open(fin):
    if l.startswith('E '): flush(); cur = {'n': int(l.split()[1])}
    elif l.startswith('A 0 GenPdfInfo'):
        f = l.split(); cur['id1'], cur['id2'], cur['x1'], cur['x2'] = int(f[3]), int(f[4]), float(f[5]), float(f[6])
flush()
np.savez(fout, evnum=np.array(ev), doublet_n0=np.array(w)); print(f'doublet_n0: <w> = {np.mean(w):.3f} (n = {len(ev)})')
