"""Independent LO normalisation check: pure-wino chargino/neutralino pair production from Pythia's own SUSY
processes (SLHA spectrum with a degenerate wino triplet, everything else decoupled), same PDF (pSet 20), 13 TeV.
usage: wino_pythia_check.py M_GeV [nev=20000]    -> prints sigma(chi+chi-), sigma(chi+chi0), sigma(chi-chi0) in fb, no acceptance cut.
Compare with the Pythia-quirk samples:  per wino  sigma_DY = sigma_gen(DY) * <w_triplet>,  sigma_W* = kappa_W^2 * sigma_gen(W*) (kappa_W^2 = 2)."""
import sys, pythia8
M = float(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
H = 5.0e4
slha = f"""BLOCK MODSEL
 1 0
BLOCK SMINPUTS
 1 1.27916e+02
 2 1.16638e-05
 3 1.18e-01
 4 9.1187e+01
 5 4.18
 6 1.73e+02
 7 1.777
BLOCK MINPAR
 3 10.
BLOCK MASS
 1000001 {H}
 1000002 {H}
 1000003 {H}
 1000004 {H}
 1000005 {H}
 1000006 {H}
 2000001 {H}
 2000002 {H}
 2000003 {H}
 2000004 {H}
 2000005 {H}
 2000006 {H}
 1000011 {H}
 1000012 {H}
 1000013 {H}
 1000014 {H}
 1000015 {H}
 1000016 {H}
 2000011 {H}
 2000013 {H}
 2000015 {H}
 1000021 {H}
 1000022 {M}
 1000023 {H}
 1000025 {H+10}
 1000035 {H+20}
 1000024 {M}
 1000037 {H+10}
 25 125.
 35 {H}
 36 {H}
 37 {H}
BLOCK NMIX
 1 1 0.
 1 2 1.
 1 3 0.
 1 4 0.
 2 1 1.
 2 2 0.
 2 3 0.
 2 4 0.
 3 1 0.
 3 2 0.
 3 3 1.
 3 4 0.
 4 1 0.
 4 2 0.
 4 3 0.
 4 4 1.
BLOCK UMIX
 1 1 1.
 1 2 0.
 2 1 0.
 2 2 1.
BLOCK VMIX
 1 1 1.
 1 2 0.
 2 1 0.
 2 2 1.
BLOCK STOPMIX
 1 1 1.
 1 2 0.
 2 1 0.
 2 2 1.
BLOCK SBOTMIX
 1 1 1.
 1 2 0.
 2 1 0.
 2 2 1.
BLOCK STAUMIX
 1 1 1.
 1 2 0.
 2 1 0.
 2 2 1.
BLOCK ALPHA
 -0.1
BLOCK HMIX Q= 1000.
 1 {H}
 2 10.
 3 2.4600e+02
 4 {H*H}
"""
open('wino.spc', 'w').write(slha)
p = pythia8.Pythia("", False)
for s in ["Beams:eCM = 13000", "PDF:pSet = 20", "SUSY:all = off", "SUSY:qqbar2chi+chi- = on", "SUSY:qqbar2chi+-chi0 = on",
          "SLHA:file = wino.spc", "SLHA:minMassSM = 100.", "PartonLevel:all = off", "HadronLevel:all = off", "ProcessLevel:resonanceDecays = off",
          "Next:numberCount = 0", "Random:setSeed = on", "Random:seed = 4242"]: p.readString(s)
p.init()
for _ in range(N):
    p.next()
info = p.infoPython() if hasattr(p, 'infoPython') else p.info
res = {}
for code in info.codesHard():
    res[code] = (info.nameProc(code), info.sigmaGen(code) * 1e12, info.sigmaErr(code) * 1e12)   # mb -> fb
for c, (n, s, e) in sorted(res.items()): print(f"M={M:g} {c} {n}: {s:.4g} +- {e:.2g} fb")
