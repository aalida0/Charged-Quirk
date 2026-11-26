# Charged-Quirk

Code for a pheno study of **SU(2)-multiplet quirks whose charged component decays**, Q+ -> Q0 + X, with the infracolor string still attached (terminating non-helical tracks). Companion paper source: `Charged-Quirk-draft`.

## Layout

| Path | Content | Status |
|---|---|---|
| `ctau/` | One-loop gauge mass splitting and Q+ -> Q0 + X decay length vs Delta m (doublet, triplet, coupling-suppressed) | first version, needs the pi-pi channel and finite-M check |
| `docs/PLAN.md` | Paper 1 plan, work packages, open checks | draft |

Planned (not started): `strings/` (string propagation with a decay step), `recast/` (public CMS HSCP / disappearing-track recasts), `plots/`.

## External code

The baseline propagator is the public QuirkPheno repository (T. A. Vami, Apache-2.0, https://github.com/tvami/QuirkPheno; derived from the ATLAS Athena Quirks Geant4 extension). It is cited, not copied here. Anything derived from it must keep its NOTICE and license.

## Run

```
cd ctau && python3 ctau.py        # prints Delta m and c*tau benchmarks
python3 plot_map.py               # writes ctau_map.png and ctau_map.csv
```

Needs numpy, scipy, matplotlib.

## Caveats on the current numbers

Static-source limit with a vector-like W* coupling; pi/e/mu below 0.5 GeV, partonic above 1.5 GeV, interpolated in between (unreliable). pi-pi channel not included. Quote at most 2 significant figures.
