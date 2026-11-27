# Scans of decaying charged quirks (patched QuirkPheno, LO, generic CMS-like detector)

Two parts.

1. **Lambda x c*tau grid at M = 500 GeV** (`make_sample.sh`, `scan.sh`, `point.sh`, `analyze_scan.py`, `plot_scan.py`; data in `data/scan_triplet_M500.csv`, `data/benchmark_*_M500.csv`). Drell-Yan Q+Q- pairs (sigma = 3.1 fb in the generated process, 3806 of 4000 events with |eta| < 2.5), triplet decay table throughout the grid, Lambda = 10, 30, 100, 300, 1000, 3000 eV, c*tau = 1, 5, 20, 54, 200, 1000, 5000 mm, **2000 charged quirks at every point** (the Lambda = 3 keV, c*tau = 5 m point was re-run with 1000 events). Same 1000 generator events at every point. Columns: see `analyze_scan.py` (`f_tracker`, `f_nl_ge3`, `dke_median`, ...). `helix_scan_M500.csv` (`analyze_helix.py`, `plot_helix.py`, `helix_fraction_M500.png`): helical-compatible fraction of the tracks with >= 4 tracker hits (circle-fit RMS < 0.1 mm, unsmeared hits).

2. **Mass scan** (`massscan/`, outputs in `data/massscan/`): M = 300, 500, 750, 1000, 1500, 2000 GeV, doublet (T = 1/2, Y = 1/2, kappa = 1) and triplet (T = 1, Y = 0, kappa = 2) with the finite-M one-loop Delta m and the decay table of each case (`massscan_params.txt`: doublet Delta m 314 -> 349 MeV, c*tau 9.5 -> 6.6 mm; triplet 164 -> 167 MeV, 59.5 -> 54.2 mm), Lambda = 30, 100, 300, 1000 eV, 1000 events per point, both production modes.
   * `massscan/gen_mass.sh M nev`: Pythia samples (DY 17 17bar; W*-mediated 17 18 / 17bar 18bar via process 3021, off-shell W' at m_W with vector-like coupling) filtered to |eta| < 2.5; `generator_summary.txt` has sigma and the kept fractions.
   * `massscan/run_ms.sh`, `ana2.sh`: simulation grid and analysis (`analyze_ms.py`, `analyze_helix.py`); `summarize_ms.py` -> `massscan_summary.csv` (cross sections and yields in fb; multiply by the luminosity in fb^-1); `plot_ms.py` -> `massscan_{xsec,yield}.png`.
   * DY events are re-weighted event by event from Pythia's gamma*/Z coupling of a chiral 4th-generation lepton to the vector-like couplings of the multiplet member (`reweight_dy.py`: <w> = 1.11-1.14 doublet, 1.45-1.53 triplet); the W* rate is multiplied by kappa_W^2 = 1 (doublet), 2 (triplet).
   * Observables per point: charged-quirk decay in the tracker, terminating-track candidate (decay in the tracker after hits in >= 3 layers), neutral-pair-leaves-detector classification with the generator-level pair pT as MET proxy (`f_evt_neutralpair_*`, `f_evt_timeout` = a quirk is still alive at the 500 ns limit), helical-compatible fraction.
   * `recast_hscp.py` -> `data/massscan/hscp_recast.csv`: probability that a charged quirk reaches r = 116 / 129 / 420 cm before decaying (straight line, no string slowing) and the surviving cross section, for the predicted c*tau and for a c*tau scan.
   * Sensitivity runs at M = 500 GeV (`massscan/extra_runs.sh`, `data/sensitivity/`): (i) doublet decay table at the triplet c*tau and Delta m: terminating-track fractions agree within 4%; (ii) ITk-like layout (`massscan/itk_like.card`, radii typed from memory, NOT verified; B = 2 T): terminating-track fraction x 0.4-0.5 for Lambda <= 300 eV and x 0.08-0.36 at 1 keV, neutral-pair yield unchanged.

## Main numbers (LO, |eta| < 2.5, no k-factor)

* Production at M = 500 GeV: DY 3.3 (doublet) / 4.3 (triplet) fb; W* 18.6 / 37 fb. W* is 6-16 times DY over 300-2000 GeV.
* Triplet, decays in the tracker after >= 3 layers: 39, 4.5, 0.36, 0.046, 2.5e-3, 2.3e-4 fb at M = 300 ... 2000 GeV (Lambda = 100 eV). Doublet: 0.37 fb at 300 GeV, 3e-3 fb at 500, none above (c*tau = 7-10 mm).
* Neutral pair leaves after a decay with pT(pair) > 250 GeV: triplet 2.5 fb (500 GeV), 0.083 fb (1 TeV), 1e-3 fb (2 TeV); independent of Lambda.
* HSCP: P(charged quirk reaches the end of the tracker) = 7.6e-5 ... 2.6e-8 (triplet, M = 300 ... 2000 GeV), < 1e-15 doublet; 0.24 at c*tau = 1 m (M = 500).
* Helical-compatible fraction (>= 4 hits): 1.00 for Lambda <= 30 eV; at c*tau = 1 (5) m 0.92 (0.76) at 100 eV, 0.49 (0.21) at 300 eV, 0.27 (0.17) at 1 keV.

## Caveats

* Leading order; the W* and DY samples are Pythia templates (process 3021 with a vector-like W'; gamma*/Z with chiral 4th-generation lepton couplings, re-weighted analytically with a massless-quark LO formula whose angular mass terms are approximate). No k-factors, no electroweak corrections, Pythia internal PDF set 20 instead of LHAPDF NNPDF3.1.
* Generic detectors only: CMS-like (the card is a reconstruction of public dimensions, see `cards/README.md` upstream) and an ITk-like variant with unverified radii. Hits are sensitive-volume visits, not reconstructed tracks; there is no smearing and no pattern recognition. The helical fraction is a proxy.
* The mass-scan yields are before any reconstruction or trigger efficiency; the neutral-pair numbers are generator-level (no ISR-jet, trigger or lepton selection).
* The per-event agreement of the handover test is limited by step discretisation (see `../quirkdecay/README.md`).
* The patched transport kills a quirk whose forced step ends outside the world; before the fix, Lambda = 10 and 30 eV, c*tau = 1 m runs segfaulted at event 530 of the 500 GeV DY sample. The underlying tiny hIoni step is not identified.
* Not done: conversion of the recasts into mass limits (needs the published limit curves and efficiency tables of CMS-EXO-18-002 and the ATLAS/CMS disappearing-track searches, which could not be fetched here); two-loop splitting; reconstruction-level terminating-track efficiency.
