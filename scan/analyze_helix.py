"""Helix-compatibility of charged-quirk tracks in the tracker (proxy for 'would a standard helical track fit accept it').
usage: analyze_helix.py <dir with L<Lambda>_c<ctau>.root> <out.csv> [cut_mm=0.1] [first_nev]
For each charged quirk (|q_pdg|=17) with >= 4 tracker hits before it decays or leaves, fit a circle to the transverse
hit mid-points (Kasa fit) and take the RMS residual.  'helical-compatible' = RMS < cut (default 100 um, ~5x the pixel
resolution).  Hit positions are the exact sensitive-volume entry/exit mid-points: no smearing, no pattern recognition."""
import sys, glob, re, csv, numpy as np, uproot
d, out = sys.argv[1], sys.argv[2]
CUT = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
NEV = int(sys.argv[4]) if len(sys.argv) > 4 else None
TRK = [2, 4, 5, 6, 7, 8]
def circ_rms(x, y):
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]; b = x ** 2 + y ** 2
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    R = np.sqrt(max(c[2] + c[0] ** 2 + c[1] ** 2, 0.))
    return np.sqrt(np.mean((np.hypot(x - c[0], y - c[1]) - R) ** 2))
rows = []
for f in sorted(glob.glob(f'{d}/L*_c*.root')):
    m = re.search(r'L([\d.]+)_c([\d.]+)\.root', f); L, ct = float(m.group(1)), float(m.group(2))
    a = uproot.open(f)['events'].arrays(['q_pdg', 'h_q', 'h_det', 'h_x0', 'h_y0', 'h_x1', 'h_y1', 'h_t0'], library='np', entry_stop=NEV)
    rms = []
    for i in range(len(a['q_pdg'])):
        for s in (0, 1):
            if abs(int(a['q_pdg'][i][s])) != 17: continue
            sel = (a['h_q'][i] == s) & np.isin(a['h_det'][i], TRK)
            if sel.sum() < 4: continue
            o = np.argsort(a['h_t0'][i][sel])
            x = ((a['h_x0'][i] + a['h_x1'][i]) / 2)[sel][o]; y = ((a['h_y0'][i] + a['h_y1'][i]) / 2)[sel][o]
            rms.append(circ_rms(x, y))
    r = np.array(rms)
    rows.append(dict(Lambda_eV=L, ctau_mm=ct, n_ge4hits=len(r), f_helical=float(np.mean(r < CUT)) if len(r) else np.nan,
                     rms_median_mm=float(np.median(r)) if len(r) else np.nan, rms_p90_mm=float(np.percentile(r, 90)) if len(r) else np.nan))
    print(f"{d.split('/')[-1]} L={L:g} ct={ct:g} n={len(r)} f_helical={rows[-1]['f_helical']:.3f} rms_med={rows[-1]['rms_median_mm']:.3g} mm", flush=True)
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
