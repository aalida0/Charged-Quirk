# Paper 1 plan: decaying charged quirks

**Claim.** For an SU(2) doublet/triplet quirk the charged component decays at tracker scale (c*tau ~ 0.6 cm doublet, ~5 cm triplet, pure EW). With an unbroken infracolor string the charged leg's non-helical track terminates while the neutral leg stays attached. No quirk paper found in INSPIRE treats this in the macroscopic-string regime (microscopic regime: Harnik et al. 1106.2569; Forsyth-Low-Tenney 2504.02940).

## Work packages

1. **Decay map (done).** Delta m and c*tau with pi-pi and finite-M. TODO: compare Forsyth Appendix A.1 normalization (their tau ~ 1e-12 s at 10 GeV is ~4e5 longer than full strength; Harnik Eq. 17 agrees with full strength).
2. **String + decay propagator (done, see quirkdecay/).** Extend the QuirkPheno transport (quirkcore) with a neutral quirk species; carry the string state across the particle swap (string state is keyed by PDG code there).
3. **Regime map (first scan done, scan/; M = 500 GeV only; mass dependence, W*-mediated production and analytic boundaries open).** Lambda vs c*tau vs oscillation/turnaround scales; adapt Harnik's decay-vs-energy-loss survival probability P = 1 - exp(-t_loss/t).
4. **Recasts (public only).** CMS HSCP (EXO-18-002, arXiv:2410.09164, HEPData ins2840007) and public disappearing-track searches, Lambda below ~100 eV.
5. **Turnaround-window selection.** Terminating coplanar-hit pattern, ~100 eV-10 keV, where the HSCP recast has no sensitivity.

## Out of scope (Paper 2/3)

Visible decays with GeV-scale splittings and mixing; forward detectors; stopped quirks; HL-LHC projections; ML reconstruction.

## Prior-art checks still open

Full texts not read: LHCb VELO quirks (2601.09023), Li et al. 1911.02223 and 2002.07503, the FASER search body. Re-run the INSPIRE audit before submission (abstract-field queries were unreliable; use free-text queries).


## Result that changes the framing (scan/)

For the pure-EW benchmarks (doublet c*tau = 6.2 mm, triplet 54 mm at M = 500 GeV) the charged quirk almost never leaves hits in three tracker layers (0.05% / 7-9%), so a reconstructed terminating non-helical track needs c*tau >~ 0.2 m. For the pure-EW benchmarks the tracker sees a quirk pair whose charged legs vanish in or before the pixel detector. The track-based signature therefore needs a longer c*tau (coupling-suppressed or mixed multiplets, Paper 2) or a different observable (neutral-pair string signature, ISR recoil). Decide this before spending effort on the disappearing-track recast.
