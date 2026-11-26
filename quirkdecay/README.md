# quirkdecay: charged-quirk decay step for QuirkPheno

A patch against [tvami/QuirkPheno](https://github.com/tvami/QuirkPheno) at commit `d91cb27` (Apache-2.0; derived from the ATLAS Athena Quirks extension; keep upstream LICENSE/NOTICE) that lets a charged quirk Q+- decay to a neutral quirk Q0 plus soft visible products while the infracolor string stays attached. Cite the QuirkPheno code and paper for everything the patch does not change.

## What it does

* `QuirkDecay` (new `G4VDiscreteProcess`): mean free path beta*gamma*c*tau, which is memoryless, so re-sampling when a suspended track resumes is exact. At decay the parent is killed, a neutral quirk (PDG |pdg|+1, or the same PDG in identity mode) and the products are created, and the string state is handed over to the neutral track (`BeginHandover()` in `QuirkInfracolorForce`), so no partner is freed.
* Kinematics: own exact 2/3-body generator (Geant4's `G4PhaseSpaceDecayChannel` fails numerically for M0 ~ 500 GeV, Delta ~ 0.36 GeV). Weights: constant for one pion; spin-averaged heavy-fermion tensor x lepton tensor for l nu (vector-like W* coupling); (k1-k2)^m (k1-k2)^n for pi pi0. No static-source shortcut: the heavy recoil is exact.
* Decay table (`make_decay_table.py <dm_GeV> <kappa> <M_GeV> <out.txt> [--charge -1]`): partial widths from `../ctau/ctau.py`; lines "weight daughter(s)", at most two visible daughters, charge conservation checked; antiquirk channels use conjugates.
* Macro commands: `/quirk/decayCtau`, `/quirk/decayDelta`, `/quirk/decayTable`, `/quirk/decayIdentity`. New ntuple columns: decay vertex/time, product kinematics and product hits (neutral quirk makes no hits).

## Build and run

```
QP_WORK=/path/to/work ./apply_and_build.sh     # clone, pin, git apply, build
cd $QP_WORK/QuirkPheno
build/quirksim/quirksim -m macros/ident.mac -t 2 -o ident.root   # identity test
build/quirksim/quirksim -m macros/lep.mac   -t 2 -s 21 -o lep.root
python3 /path/to/quirkdecay/validate_decay.py lep.root 500 0.356 6.21196
```

`apply_and_build.sh` was not run end-to-end from a fresh clone; the patch itself was checked with `git apply --check` on a pristine checkout of `d91cb27`, and the build was done from the equivalent working tree with Geant4 11.4.3.

## Validation done (M = 500 GeV, Delta = 356 MeV, c*tau = 6.212 mm, leptonic table, 2 x 500 events)

* Proper decay time vs exponential: KS p = 0.04 and 0.40 (mean 6.27 +- 0.21 and 6.39 +- 0.20 mm).
* Four-momentum conserved per decay to float precision.
* Rest-frame mean energy: e 173.2/181.0 (+-3.1) vs 177.9 MeV; mu 207.8/208.5 (+-2.4) vs 206.8 MeV (static-source spectrum E p (Delta-E)^2).
* Identity mode (neutral daughter same PDG): string handover reproduces the undecayed run (see `macros/ident.mac`).
* An earlier flat-invariant-matrix-element generator gave a 3.7% high muon mean (214 vs 207 MeV); the exact flat-Dalitz MC agreed with it (216 +- 1), showing the static reference needs the E_l E_nu matrix-element factor. The tensor weighting fixes this.
* pi pi0 weight: non-negative and below its rejection bound on random configurations; the pi pi0 channel has not been validated against its analytic spectrum in the simulation.

## Limits

* Drell-Yan Q+Q- only: no W-mediated or neutral-pair production in the Pythia card; the neutral PDG is accepted in gun/HepMC mode only.
* Hadronic continuum above ~0.8 GeV is not modelled; at most two visible products; rho dynamics absent from the pi pi0 weight.
* Stopped quirks are killed before they decay unless `keepStopped` is set.
* The neutral quirk has string interactions only.
* The handover is validated by the identity test only, not against external code.
* Delta m(M) is the one-loop Dirac-doublet/triplet value; no two-loop correction (a few MeV).
