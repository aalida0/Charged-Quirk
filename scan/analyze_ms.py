"""Mass-scan analysis of quirksim runs with the quirkdecay patch (DY Q+Q- or W*-mediated Q+Q0 samples).
usage: analyze_ms.py <dir with L<Lambda>_c<ctau>.root> <out.csv> <M_GeV> [met_cut_GeV=250]
Per CHARGED quirk (|q_pdg| = 17; the neutral partner of a W* pair is excluded): decay inside the tracker, layers hit
before the decay, kinetic-energy loss.  Per event: terminating-track candidate (a charged quirk decays inside the tracker
after hits in >= 3 layers), and the neutral-pair signature (the pair leaves the detector, i.e. end reason 1 or 2 for the first
quirk and 1, 2 or 4 (ended because its partner did) for the other, after a charged-quirk decay; 'timeout' = still alive at the 500 ns limit), with the generator-level pair pT as the MET proxy (the decay products are soft)."""
import sys, glob, re, csv, numpy as np, uproot
d, out, M = sys.argv[1], sys.argv[2], float(sys.argv[3]) * 1e3
METCUT = float(sys.argv[4]) * 1e3 if len(sys.argv) > 4 else 250e3
TRK = {2, 4, 5, 6, 7, 8}
R_TRK, Z_TRK, R_PIX = 1160., 2800., 160.
ESC, PARTNER, TIMEOUT = {1, 2}, 4, 6   # end reasons: envelope/world exit; ended because the partner ended; 500 ns time limit
rows, evdump = [], {}
for f in sorted(glob.glob(f'{d}/L*_c*.root')):
    m = re.search(r'L([\d.]+)_c([\d.]+)\.root', f); L, ct = float(m.group(1)), float(m.group(2))
    a = uproot.open(f)['events'].arrays(['dec_side', 'dec_x', 'dec_y', 'dec_z', 'dec_pe', 'q_pdg', 'q_e0', 'q_x0', 'q_y0', 'q_z0',
        'q_px0', 'q_py0', 'q_reason', 'h_q', 'h_det', 'h_layer', 'genEvt'], library='np')
    n = len(a['q_pdg']); rec = []
    flags = []
    ev_term = ev_neutral = ev_neutral_met = ev_ann = ev_esc = ev_dec = 0
    for i in range(n):
        side_dec = {int(s): j for j, s in enumerate(a['dec_side'][i])}
        hq, hd, hl = a['h_q'][i], a['h_det'][i], a['h_layer'][i]
        any_term = False
        for s in (0, 1):
            if abs(int(a['q_pdg'][i][s])) != 17: continue
            ke0 = a['q_e0'][i][s] - M
            if s in side_dec:
                j = side_dec[s]; x, y, z = a['dec_x'][i][j], a['dec_y'][i][j], a['dec_z'][i][j]; ke = a['dec_pe'][i][j] - M; dec = 1
            else:
                x = y = z = np.nan; ke = np.nan; dec = 0
            sel = (hq == s) & np.isin(hd, list(TRK))
            nl = len({(int(dd), int(ll)) for dd, ll in zip(hd[sel], hl[sel])})
            r = np.hypot(x, y) if dec else np.nan
            Ld = np.sqrt((x - a['q_x0'][i][s])**2 + (y - a['q_y0'][i][s])**2 + (z - a['q_z0'][i][s])**2) if dec else np.nan
            intrk = bool(dec and r < R_TRK and abs(z) < Z_TRK)
            rec.append((dec, intrk, bool(dec and r < R_PIX), nl, (ke0 - ke) / ke0 if dec else np.nan, Ld, ke0, ke0 - ke if dec else np.nan))
            if intrk and nl >= 3: any_term = True
        rs = [int(v) for v in a['q_reason'][i]]
        ndec = len(side_dec)
        pt = np.hypot(a['q_px0'][i].sum(), a['q_py0'][i].sum())
        esc = all(v in ESC or v == PARTNER for v in rs) and any(v in ESC for v in rs); ann = any(v == TIMEOUT for v in rs)
        ev_term += any_term; ev_dec += ndec > 0; ev_esc += esc; ev_ann += ann
        ev_neutral += (ndec > 0 and esc); ev_neutral_met += (ndec > 0 and esc and pt > METCUT)
        flags.append((int(a['genEvt'][i]), any_term, ndec > 0 and esc, ndec > 0 and esc and pt > METCUT, any(v == TIMEOUT for v in rs), pt))
    evdump[f'L{L:g}'] = np.array(flags, float)   # columns: genEvt, termtrack, neutralpair_escape, neutralpair_met, timeout, pair pT (MeV)
    r = np.array(rec, float); dec = r[:, 0] == 1; nq = len(r)
    rows.append(dict(Lambda_eV=L, ctau_mm=ct, nevents=n, nquirks_charged=nq, f_decayed=dec.mean(), f_tracker=r[:, 1].mean(),
        f_pixel=r[:, 2].mean(), f_nl_ge3=np.mean(dec & (r[:, 3] >= 3)), f_nl_ge5=np.mean(dec & (r[:, 3] >= 5)), f_survive=1 - dec.mean(),
        dke_median=np.nanmedian(r[dec, 4]) if dec.any() else np.nan, dE_med_GeV=np.nanmedian(r[dec, 7]) / 1e3 if dec.any() else np.nan,
        Ld_med_mm=np.nanmedian(r[dec, 5]) if dec.any() else np.nan,
        f_evt_termtrack=ev_term / n, f_evt_anydecay=ev_dec / n, f_evt_escape=ev_esc / n, f_evt_timeout=ev_ann / n,
        f_evt_neutralpair_escape=ev_neutral / n, f_evt_neutralpair_met=ev_neutral_met / n))
    print(f"{d.split('/')[-1]} L={L:g} ct={ct:g} n={n} termtrack={ev_term/n:.3f} neutralpair_escape={ev_neutral/n:.3f} timeout={ev_ann/n:.3f}", flush=True)
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
np.savez(out.replace('.csv', '_ev.npz'), **evdump)
