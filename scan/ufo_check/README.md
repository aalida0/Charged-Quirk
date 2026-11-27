# Cross-check of the QuirkSim Dirac-triplet NLO-QCD UFO model (DiracTripletEW_NLO_UFO)

Model provided by the user (provenance file: loop_sm base, NLO QCD only, no hidden-colour index, no mass splitting, no decay). Run with MG5_aMC 3.8.1, 13 TeV, restriction `-no_b_mass`.

* `tri_lo.mg5`, `run_tri_lo.sh`: LO, six channels (q1+q1-, q2-q2+, q1+q0~, q1-q0, q2-q0~, q2+q0), nn23lo1.
* `tri_nlo.mg5`, `run_tri_nlo.sh`: fixed-order NLO QCD of q1+q1-, q1+q0~, q1-q0 (= the wino-equivalent total; the Dirac triplet is twice this), nn23nlo, CutTools+IREGI (no ninja/collier: no network).
* Result table: `../data/massscan/ufo_nlo_crosscheck.csv`.

Findings
1. Matrix elements: with identical EW inputs (aEWM1=127.934, same Gmu-derived MW) and a fixed scale the u ubar -> pair cross section equals that of the MSSM_SLHA2 pure wino to 0.1% (M=1000: 1.501e-4 pb both). The 6-10% difference at the default dynamical scale is the default scale choice, not the model. The Dirac triplet is exactly twice the wino per channel.
2. NLO QCD (half-triplet = wino-equivalent total) over the public NLO+NLL wino tables: 0.970 (300), 0.970 (500), 0.952 (1000 GeV) at aEWM1=127.934; 0.997 at 500 GeV with the model default aEWM1=132.507 (Gmu scheme). Scale envelope +-2-3.6%. The NLL resummation is expected to add a few per cent at these masses, so the model reproduces the public normalisation within the EW-input scheme uncertainty (~3%).
3. Notes: object_library.py imports `six`, which MG5 cannot find because it replaces sys.path during the model import (patched in the local sandbox copy only: `six.iteritems(d)` -> `d.items()`; original untouched). The banner lists `9000113 2.25e+02` for q2- but the LHE masses and cross sections use MQC (checked). The transport (quirksim PrimaryGenerator::isQuirkPdg) identifies quirks by one PDG q and q+1; the UFO PDGs 9000111/9000112/9000113 need an extended list.
