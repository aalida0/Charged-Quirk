# Hard-process-only sigma for the quirk DY card vs Pythia internal PDF set. usage: pdf_xsec_check.py card.cmnd <mass> [PDF:pSet=NN ...]
# (card = QuirkPheno card with the HadronLevel:QED and PDF:pSet lines removed; valid for the taup/doublet cards at m<=~600; fails at higher m because of the Z-prime resonance-decay setup)
import sys, pythia8
card, m = sys.argv[1], float(sys.argv[2]); extra = sys.argv[3:]
p = pythia8.Pythia("", False)
assert p.readFile(card)
for s in extra + ["17:m0=%f"%m,"17:mWidth=0.","17:mayDecay=off","17:isResonance=false","32:mMin=%f"%(2*m+1),
          "PartonLevel:all=off","HadronLevel:all=off","ProcessLevel:resonanceDecays=on","Next:numberCount=0","Random:setSeed=on","Random:seed=3"]:
    p.readString(s)
p.init()
for _ in range(30000): p.next()
print("%.5e +- %.1e pb"%(p.infoPython().sigmaGen()*1e9, p.infoPython().sigmaErr()*1e9))
