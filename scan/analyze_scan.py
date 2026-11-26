"""Summarise a Lambda x c*tau scan of quirksim runs (quirkdecay patch, DY pairs, triplet decay table).
usage: analyze_scan.py <dir with L<Lambda>_c<ctau>.root> <out.csv> [M_GeV]
Decay distance is measured from the quirk's production point (beam-spot sigma_z ~ 35 mm in the generated sample).
Per charged quirk (both signs of the pair): decay inside the tracker, layers hit before the decay, kinetic-energy
loss to the string/ionisation at the decay (string-induced non-helical bending is tied to it)."""
import sys, glob, re, numpy as np, uproot
d, out = sys.argv[1], sys.argv[2]
M = float(sys.argv[3]) * 1e3 if len(sys.argv) > 3 else 5e5       # MeV
TRK = {2, 4, 5, 6, 7, 8}                                          # PXB PXF TIB TOB TID TEC (det ids of the card)
R_TRK, Z_TRK, R_PIX = 1160., 2800., 160.                          # mm
rows = []
for f in sorted(glob.glob(f'{d}/L*_c*.root')):
    m = re.search(r'L([\d.]+)_c([\d.]+)\.root', f); L, ct = float(m.group(1)), float(m.group(2))
    a = uproot.open(f)['events'].arrays(['dec_side','dec_x','dec_y','dec_z','dec_pe','q_e0','q_x0','q_y0','q_z0',
        'h_q','h_det','h_layer','q_dec'], library='np')
    n = len(a['q_dec'])
    rec = []
    for i in range(n):
        side_dec = {int(s): j for j, s in enumerate(a['dec_side'][i])}
        hq, hd, hl = a['h_q'][i], a['h_det'][i], a['h_layer'][i]
        for s in (0, 1):
            e0 = a['q_e0'][i][s]; ke0 = e0 - M
            if s in side_dec:
                j = side_dec[s]
                x, y, z = a['dec_x'][i][j], a['dec_y'][i][j], a['dec_z'][i][j]
                ke = a['dec_pe'][i][j] - M
                dec = 1
            else:
                x = y = z = np.nan; ke = np.nan; dec = 0
            sel = (hq == s) & np.isin(hd, list(TRK))
            nl = len({(int(dd), int(ll)) for dd, ll in zip(hd[sel], hl[sel])})
            rec.append((dec, np.hypot(x, y), abs(z) if dec else np.nan, (ke0 - ke) / ke0 if dec else np.nan, nl,
                        np.sqrt((x - a['q_x0'][i][s])**2 + (y - a['q_y0'][i][s])**2 + (z - a['q_z0'][i][s])**2) if dec else np.nan, ke0, (ke0 - ke) if dec else np.nan))
    r = np.array(rec, float); dec = r[:, 0] == 1
    intrk = dec & (r[:, 1] < R_TRK) & (r[:, 2] < Z_TRK)
    nq = len(r)
    rows.append(dict(Lambda_eV=L, ctau_mm=ct, nquirks=nq, f_decayed=dec.mean(), f_tracker=intrk.mean(),
        f_pixel=(dec & (r[:, 1] < R_PIX)).mean(), nl_mean_decayed=np.nanmean(r[dec, 4]) if dec.any() else np.nan,
        f_nl_ge3=np.mean(dec & (r[:, 4] >= 3)), f_nl_ge5=np.mean(dec & (r[:, 4] >= 5)),
        f_survive=1 - dec.mean(), dke_median=np.nanmedian(r[dec, 3]) if dec.any() else np.nan,
        dke_p90=np.nanpercentile(r[dec, 3], 90) if dec.any() else np.nan,
        f_dke_gt5=np.mean(dec & (r[:, 3] > 0.05)),
        Ld_med_mm=np.nanmedian(r[dec, 5]) if dec.any() else np.nan, ke0_med_GeV=np.nanmedian(r[:, 6]) / 1e3,
        dE_med_GeV=np.nanmedian(r[dec, 7]) / 1e3 if dec.any() else np.nan))
    print(f"L={L:g} ct={ct:g} nq={nq} dec={dec.mean():.3f} trk={intrk.mean():.3f} nl>=3={rows[-1]['f_nl_ge3']:.3f} dke_med={rows[-1]['dke_median']:.3g}", flush=True)
import csv
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
