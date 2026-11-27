"""Smeared-hit helix-compatibility proxy.
usage: analyze_helix_smeared.py <dir with L<Lambda>_c<ctau>.root> <out.csv> [chi2ndf_cut=3] [first_nev] [seed]
Same hits as analyze_helix.py, but every hit mid-point is smeared by a Gaussian in the transverse plane
(sigma_pix = 10 um for PXB/PXF, sigma_strip = 30 um for TIB/TOB/TID/TEC, isotropic) and the circle is fitted by an
algebraic (Kasa) fit followed by a few Gauss-Newton iterations of the geometric chi2 with those sigmas.  A charged
quirk with >= 4 tracker hits is 'helix-compatible' if chi2/ndf (ndf = n_hits - 3) is below the cut, AND (to mimic the
radius/pT requirement of a standard fit) |R| > 0.5 m, i.e. pT > ~0.15 GeV at 3.8 T.  Pattern recognition, material
effects and a longitudinal (rz) fit are not modelled; the number is a proxy, not a reconstruction efficiency."""
import sys, glob, re, csv, numpy as np, uproot
d, out = sys.argv[1], sys.argv[2]
CUT = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
NEV = int(sys.argv[4]) if len(sys.argv) > 4 else None
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 12345)
TRK = [2, 4, 5, 6, 7, 8]; PIX = {2, 4}
def fit(x, y, s):
    w = 1. / s ** 2
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]; b = x ** 2 + y ** 2
    c, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], b * np.sqrt(w), rcond=None)
    xc, yc = c[0], c[1]; R = np.sqrt(max(c[2] + xc ** 2 + yc ** 2, 1e-12))
    for _ in range(5):                                   # geometric Gauss-Newton on (xc, yc, R)
        dx, dy = x - xc, y - yc; rr = np.hypot(dx, dy) + 1e-12
        res = rr - R
        J = np.c_[-dx / rr, -dy / rr, -np.ones_like(x)] / s[:, None]
        step, *_ = np.linalg.lstsq(J, -res / s, rcond=None)
        xc, yc, R = xc + step[0], yc + step[1], R + step[2]
    res = (np.hypot(x - xc, y - yc) - R) / s
    return float(np.sum(res ** 2)) / (len(x) - 3), abs(R)
rows = []
for f in sorted(glob.glob(f'{d}/L*_c*.root')):
    m = re.search(r'L([\d.]+)_c([\d.]+)\.root', f); L, ct = float(m.group(1)), float(m.group(2))
    a = uproot.open(f)['events'].arrays(['q_pdg', 'h_q', 'h_det', 'h_x0', 'h_y0', 'h_x1', 'h_y1', 'h_t0'], library='np', entry_stop=NEV)
    c2, Rs = [], []
    for i in range(len(a['q_pdg'])):
        for s in (0, 1):
            if abs(int(a['q_pdg'][i][s])) != 17: continue
            sel = (a['h_q'][i] == s) & np.isin(a['h_det'][i], TRK)
            if sel.sum() < 4: continue
            o = np.argsort(a['h_t0'][i][sel])
            x = ((a['h_x0'][i] + a['h_x1'][i]) / 2)[sel][o]; y = ((a['h_y0'][i] + a['h_y1'][i]) / 2)[sel][o]
            sg = np.where(np.isin(a['h_det'][i][sel][o], list(PIX)), 0.010, 0.030)
            x = x + rng.normal(0, sg); y = y + rng.normal(0, sg)
            chi, R = fit(x, y, sg); c2.append(chi); Rs.append(R)
    c2, Rs = np.array(c2), np.array(Rs)
    ok = (c2 < CUT) & (Rs > 500.) if len(c2) else np.array([], bool)
    rows.append(dict(Lambda_eV=L, ctau_mm=ct, n_ge4hits=len(c2), f_helical_smeared=float(np.mean(ok)) if len(c2) else np.nan,
                     chi2ndf_median=float(np.median(c2)) if len(c2) else np.nan))
    print(f"{d.split('/')[-1]} L={L:g} ct={ct:g} n={len(c2)} f_hel_smeared={rows[-1]['f_helical_smeared']:.3f} chi2/ndf_med={rows[-1]['chi2ndf_median']:.3g}", flush=True)
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
