# Cross-check of the QuirkSim Dirac-triplet NLO-QCD UFO model (DiracTripletEW_NLO_UFO)

Model provided by the user (provenance file: loop_sm base, NLO QCD only, no hidden-colour index, no mass splitting, no decay). Run with MG5_aMC 3.8.1, 13 TeV, restriction `-no_b_mass`.

* `tri_lo.mg5`, `run_tri_lo.sh`: LO, six channels (q1+q1-, q2-q2+, q1+q0~, q1-q0, q2-q0~, q2+q0), nn23lo1.
* `tri_nlo.mg5`, `run_tri_nlo.sh`: fixed-order NLO QCD of q1+q1-, q1+q0~, q1-q0 (= the wino-equivalent total; the Dirac triplet is twice this), nn23nlo, CutTools+IREGI (no ninja/collier: no network).
* Result table: `../data/massscan/ufo_nlo_crosscheck.csv`.

Findings
1. Matrix elements: with identical EW inputs (aEWM1=127.934, same Gmu-derived MW) and a fixed scale the u ubar -> pair cross section equals that of the MSSM_SLHA2 pure wino to 0.1% (M=1000: 1.501e-4 pb both). The 6-10% difference at the default dynamical scale is the default scale choice, not the model. The Dirac triplet is exactly twice the wino per channel.
2. NLO QCD (half-triplet = wino-equivalent total) over the public NLO+NLL wino tables: 0.970 (300), 0.970 (500), 0.952 (1000 GeV) at aEWM1=127.934; 0.997 at 500 GeV with the model default aEWM1=132.507 (Gmu scheme). Scale envelope +-2-3.6%. The NLL resummation is expected to add a few per cent at these masses, so the model reproduces the public normalisation within the EW-input scheme uncertainty (~3%).
3. Notes: object_library.py imports `six`, which MG5 cannot find because it replaces sys.path during the model import (patched in the local sandbox copy only: `six.iteritems(d)` -> `d.items()`; original untouched). The banner lists `9000113 2.25e+02` for q2- but the LHE masses and cross sections use MQC (checked). The transport (quirksim PrimaryGenerator::isQuirkPdg) identifies quirks by one PDG q and q+1; the UFO PDGs 9000111/9000112/9000113 need an extended list.

## Closure of the UFO sample through the existing chain (M=500 GeV, triplet)
No change of the transport: `ufo_lhe_to_quirk.py` maps the UFO ids to the ids the transport and the Pythia chain already use (charge -1 -> +17, charge +1 -> -17; the neutral partner of a W*-mediated pair takes the opposite sign, as for the W' -> -17 18 of the template chain; Q1/Q2 species dropped, degenerate). Chain: MG5 LO LHE (20000 events, all six channels, sigma = 0.0941 pb at aEWM1=127.934) -> `ufo_lhe_to_quirk.py` -> `ufo_rec_gen.sh` (quirkgen with `ufo.cmnd` = the DY Pythia card with `Beams:frameType=4`: shower, MPI, hadronisation as before; pair pT>180 GeV, |eta|<2.5 filter) -> `ufo_rec_sim.sh` (Geant4 transport, same decay table and cms_like card). Class fractions in the LHE: DY 32%, W+ 48%, W- 20% (analytic 32/48/20%).

* Recoil: fraction passing the filter, UFO+Pythia vs the Pythia template: DY 10.0+-0.4% vs 11.0+-0.2% (-8%, 2 sigma), W* 10.4+-0.3% vs 10.6+-0.2%; pair-pT median 272/270 vs 274/272 GeV, KS p = 0.68 (DY), 0.27 (W*). The template DY is a 4th-generation-lepton gamma*/Z process re-weighted analytically; the UFO is the reference.
* Transport (`closure_transport.py`, `closure_transport_M500.txt`): decayed charged endpoints per event 2.000/1.000 (DY/W*), fractions with >=3 distinct layers 0.102+-0.008 (DY) and 0.103+-0.008 (W*) at 100 eV against 0.100+-0.005 and 0.103+-0.007 for the template, and 0.018+-0.004, 0.016+-0.003 against 0.019+-0.002, 0.020+-0.003 at 1 keV. The PDG handling is therefore validated and no change of `isQuirkPdg` is needed.
* Not done here: matrix-element-merged 0/1/2-jet ISR (MG5 + Pythia) with the one-jet check at the UFO level, NLO+PS, and the scaling of the closure to the other masses.
