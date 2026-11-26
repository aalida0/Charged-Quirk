# Lambda x c*tau scan of decaying charged quirks

Chain: `make_sample.sh` (Pythia Drell-Yan Q+Q-, M = 500 GeV, |eta| < 2.5, 3806 of 4000 events kept, sigma = 3.0 fb in the generated process) -> `scan.sh "<Lambda eV list>" "<ctau mm list>" <nev>` (runs `point.sh` with the patched QuirkPheno, CMS-like card, triplet decay table, 1000 events per point = 2000 charged quirks) -> `analyze_scan.py` -> `plot_scan.py`. Paths in the shell scripts are those of the sandbox build.

Results: `data/scan_triplet_M500.csv` (Lambda = 10, 30, 100, 300, 1000, 3000 eV; c*tau = 1, 5, 20, 54, 200, 1000, 5000 mm; the Lambda = 3 keV, c*tau = 5 m point has 400 events), `data/benchmark_{doublet,triplet}_M500.csv` (c*tau = 6.21 mm with Delta m = 356 MeV, and 53.9 mm with Delta m = 167 MeV, kappa = 2).

Same 1000 generator events are used at every point, so differences between points are not statistical fluctuations of the sample. The binomial error on a fraction f from 2000 quirks is sqrt(f(1-f)/2000) <= 1.1%; the two quirks of an event are not independent.

## What the columns mean

`f_tracker`: charged quirk decays inside r < 116 cm, |z| < 280 cm. `f_pixel`: r < 16 cm. `f_nl_ge3`, `f_nl_ge5`: it decays after leaving sensitive hits in at least 3 (5) distinct tracker layers (barrel, endcap disks counted as layers). `dke_median`: median fractional loss of kinetic energy from production to the decay. `Ld_med_mm`, `dE_med_GeV`: median decay distance from the production point and median kinetic-energy loss. `f_survive`: no decay before leaving the world.

## Caveats

* Production is leading-order Drell-Yan only; no W-mediated Q+Q0 production, no EW corrections.
* `ctau` is an input, not computed per point: the triplet table (97.5% pi, 2.3% e, 0.3% mu) is used throughout the grid; only the benchmark points use the doublet table.
* Hits are the sensitive-volume visits recorded by quirksim (one per volume), not reconstructed tracks; "layers hit" is not a track-reconstruction efficiency.
* The patched transport kills a quirk whose forced step ends outside the world (see `../quirkdecay/README.md`); before the fix, the Lambda = 10 and 30 eV, c*tau = 1 m runs segfaulted at event 530 of this sample.
