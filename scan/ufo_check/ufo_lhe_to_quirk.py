#!/usr/bin/env python3
"""Remap the PDG ids of the DiracTripletEW_NLO_UFO quirks in an LHE file to the ids of the transport (|pdg| 17 charged, 18 neutral partner).
Charge convention of the transport/Pythia: +17 has charge -1 (tau' like), so charge -1 -> +17, charge +1 -> -17; the neutral partner of a W*-mediated
pair gets the opposite sign of its charged partner (as the W' -> -17 18 of the Pythia chain). Species (Q1 vs Q2) is dropped (degenerate, same transport).
UFO ids: 9000111 q1+ (charge +1), 9000112 q0, 9000113 q2- (charge -1); anti-particles carry the opposite sign.
usage: ufo_lhe_to_quirk.py in.lhe[.gz] out.lhe"""
import sys, gzip, re
fin, fout = sys.argv[1:3]
CH = {9000111: +1, -9000111: -1, 9000113: -1, -9000113: +1}   # charge of the (anti)particle
NEU = {9000112, -9000112}
opn = gzip.open if fin.endswith('.gz') else open
n = nev = 0
with opn(fin, 'rt') as f, open(fout, 'w') as o:
    ev = None
    for line in f:
        if line.startswith('<event>'): ev = [line]; continue
        if ev is None: o.write(line); continue
        ev.append(line)
        if line.startswith('</event>'):
            body = ev[2:-1]; rows = []
            for k, l in enumerate(body):
                t = l.split()
                if len(t) >= 13 and not l.startswith('<') and not l.startswith('#'):
                    rows.append((k, int(t[0]), int(t[1])))
            newid = {}
            chg = [(k, p) for k, p, st in rows if st == 1 and p in CH]
            for k, p in chg: newid[k] = -17 if CH[p] > 0 else 17            # charge -1 -> +17
            sign_ch = [newid[k] for k, _ in chg]
            for k, p, st in rows:
                if st == 1 and p in NEU:
                    assert len(sign_ch) == 1, 'neutral without exactly one charged partner'
                    newid[k] = 18 if sign_ch[0] < 0 else -18                  # opposite sign to the charged partner
            for k, l in enumerate(body):
                if k in newid:
                    t = l.split(None, 1); body[k] = f' {newid[k]:>9d} ' + t[1] if False else re.sub(r'^\s*-?\d+', f'{newid[k]:>10d}', l, count=1)
            o.write(ev[0] + ev[1]); o.writelines(body); o.write(ev[-1]); ev = None; nev += 1
print(f'remapped {nev} events')
