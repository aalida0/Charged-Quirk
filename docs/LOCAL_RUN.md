# Running QuirkPheno locally (no CVMFS / lxplus)

`local/setup_env.sh` builds QuirkPheno (T. A. Vámi, Apache-2.0) against conda-forge Geant4 11.4.3 + Pythia 8.312.
`local/run_local.sh <taup|quirk> <m> <Lambda_eV> <n> <seed> <outdir>` runs quirkgen -> prep_gen.py -> quirksim.
Upstream is validated with LCG_110 (Geant4 11.4.1); the build compiles with warnings only, but no numerical
comparison against an LCG build has been done.

## Deviations from the upstream Pythia cards (applied by sed at run time; upstream files untouched)

| Upstream line | Local | Why |
|---|---|---|
| `HadronLevel:QED = on` | removed | not recognised by Pythia 8.312 (settings lookup fails, card is rejected) |
| `PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118` | `PDF:pSet = 20` (NNPDF3.1 QCD+LUXQED NNLO, alpha_s 0.118) | LHAPDF is installed but the grid download hosts (lhapdfsets.web.cern.ch, CVMFS) are refused from this sandbox; set 20 is the closest Pythia-internal set |

PDF sensitivity of the hard cross-section (DY quirk pair, m = 500 GeV, hard process only, 30k events):

| pSet | description | sigma (pb) |
|---|---|---|
| 13 | NNPDF2.3 QCD+QED LO, as 0.130 | 2.63e-3 |
| 17 | NNPDF3.1 LO, as 0.130 | 3.48e-3 |
| 18 | NNPDF3.1 LO, as 0.118 | 3.64e-3 |
| 19 | NNPDF3.1 LUXQED NLO | 3.10e-3 |
| **20** | **NNPDF3.1 LUXQED NNLO (default here)** | **3.04e-3** |

Spread ~ +-15% around set 20 for the LO sets, ~2% between 19 and 20. The remaining difference to the true
`NNPDF31_nnlo_as_0118` (QCD-only, no photon PDF) is not quantified; for absolute rates, install the grid
(`lhapdf install NNPDF31_nnlo_as_0118`) on a machine that can reach the LHAPDF mirror and use the unpatched card.
The efficiencies (ratios) are much less PDF-sensitive than absolute rates.

## What was and was not validated

Done: full chain runs for `taup` (m = 1000, Lambda = 1 eV, 400 generated, 317 kept by `prep_gen.py`) and the
doublet `quirk` card (m = 500, Lambda = 300 eV, 30 generated).

Acceptance-level sanity check (tau', m = 1000, Lambda = 1 eV): fraction of generated events with >= 1 quirk with
pT > 55 GeV, |eta| < 1 that leaves the detector is 0.61, versus the published EXO-18-002 total efficiency 0.37
(= card_yields.json total / (sigma_ref x 101 fb^-1)). Consistent as an upper bound (the published number adds the
trigger and the identification cuts); this is **not** a closure test.

Not done (blocked in this sandbox): the full emulation + limit closure (`recast/emu_all.sh`,
`recast/closure_limits.sh`). It needs the EXO-18-002 HEPData record (trigger-efficiency tables read by
`qreco.load_trigger`; hepdata.net is refused by the proxy), ROOT and Combine (CVMFS-only in upstream). To close
the loop, run `hepdata/fetch.sh` on a machine with network access and copy `hepdata/original/` here.
