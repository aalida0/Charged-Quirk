# Paper 1 plan: decaying charged quirks

**Claim.** For an SU(2) doublet/triplet quirk the charged component decays at tracker scale (c*tau ~ 0.6 cm doublet, ~5 cm triplet, pure EW). With an unbroken infracolor string the charged leg's non-helical track terminates while the neutral leg stays attached. No quirk paper found in INSPIRE treats this in the macroscopic-string regime (microscopic regime: Harnik et al. 1106.2569; Forsyth-Low-Tenney 2504.02940).

## Work packages

1. **Decay map (done, v0).** Delta m and c*tau. TODO: pi-pi channel (open above 0.28 GeV), finite-M corrections, compare Forsyth Appendix A.1 normalization (their tau ~ 1e-12 s at 10 GeV is ~4e5 longer than full strength; Harnik Eq. 17 agrees with full strength).
2. **String + decay propagator.** Extend the QuirkPheno transport (quirkcore) with a neutral quirk species; carry the string state across the particle swap (string state is keyed by PDG code there).
3. **Regime map.** Lambda vs c*tau vs oscillation/turnaround scales; adapt Harnik's decay-vs-energy-loss survival probability P = 1 - exp(-t_loss/t).
4. **Recasts (public only).** CMS HSCP (EXO-18-002, arXiv:2410.09164, HEPData ins2840007) and public disappearing-track searches, Lambda below ~100 eV.
5. **Turnaround-window selection.** Terminating coplanar-hit pattern, ~100 eV-10 keV, where the HSCP recast has no sensitivity.

## Out of scope (Paper 2/3)

Visible decays with GeV-scale splittings and mixing; forward detectors; stopped quirks; HL-LHC projections; ML reconstruction.

## Prior-art checks still open

Full texts not read: LHCb VELO quirks (2601.09023), Li et al. 1911.02223 and 2002.07503, the FASER search body. Re-run the INSPIRE audit before submission (abstract-field queries were unreliable; use free-text queries).
