"""C5 diagnostic: per decayed charged quirk, compare the string-induced curvature shift from (a) the bending-angle proxy of atlas_dt_recast.py
(direction change between production and the last hit at r<125 mm minus the magnetic expectation) and (b) a circle fit through the production
point and the true hits at r<125 mm.  Reported as the equivalent shift of q/pT at the ATLAS field (2 T), in 1/TeV, against the tracklet
resolution sigma = 14-21/TeV.   usage: hitfit_vs_proxy.py <ms dir> <M>"""
import sys, glob, numpy as np, uproot
base, M = sys.argv[1], int(sys.argv[2]); BSIM, BATL = 3.8, 2.0
for L in (30, 100, 300, 1000):
    dp, df, nh, npix = [], [], [], 0
    for p in ('dy', 'w'):
        fn = glob.glob(f'{base}/out3/triplet_{p}_M{M}/L{L}_c*.root')[0]
        a = uproot.open(fn)['events'].arrays(['q_pdg', 'q_x0', 'q_y0', 'q_px0', 'q_py0', 'dec_side', 'h_q', 'h_det', 'h_x0', 'h_y0', 'h_x1', 'h_y1', 'h_px', 'h_py', 'h_t0'], library='np')
        for i in range(len(a['q_pdg'])):
            for j, s in enumerate(a['dec_side'][i]):
                s = int(s); pdg = int(a['q_pdg'][i][s])
                if abs(pdg) != 17: continue
                sel = (a['h_q'][i] == s) & np.isin(a['h_det'][i], [2, 4, 5, 6, 7, 8])
                if sel.sum() < 2: continue
                o = np.argsort(a['h_t0'][i][sel]); x = ((a['h_x0'][i] + a['h_x1'][i]) / 2)[sel][o]; y = ((a['h_y0'][i] + a['h_y1'][i]) / 2)[sel][o]
                xs, ys = x - a['q_x0'][i][s], y - a['q_y0'][i][s]; r = np.hypot(xs, ys); m = np.where(r < 125.)[0]
                if len(m) < 2: continue
                jj = m[-1]; px = a['h_px'][i][sel][o][jj]; py = a['h_py'][i][sel][o][jj]; p0x, p0y = a['q_px0'][i][s], a['q_py0'][i][s]
                pt0 = np.hypot(p0x, p0y) / 1e3; Lm = r[jj] / 1e3; sq = 1. if pdg > 0 else -1.
                d = (np.arctan2(py, px) - np.arctan2(p0y, p0x) + np.pi) % (2 * np.pi) - np.pi
                dphis = d - sq * 0.29979 * BSIM * Lm / pt0
                A = np.c_[xs[m], ys[m]]; b = 0.5 * (xs[m] ** 2 + ys[m] ** 2); c = np.linalg.lstsq(A, b, rcond=None)[0]; R2 = c @ c; ph0 = np.arctan2(p0y, p0x)
                kap = (-c[0] * np.sin(ph0) + c[1] * np.cos(ph0)) / R2 * 1e3
                fit = (kap - sq * 0.29979 * BSIM / pt0) * Lm
                dp.append(1e3 * dphis / (0.29979 * BATL * Lm)); df.append(1e3 * fit / (0.29979 * BATL * Lm)); nh.append(len(m))
    dp, df = np.array(dp), np.array(df); good = np.abs(dp) < 1e4
    print(f'M={M} L={L:5d} eV  N={len(dp)} <nhits>={np.mean(nh):.2f}  rms(proxy)={np.sqrt(np.mean(dp**2)):9.2f}  rms(fit)={np.sqrt(np.mean(df**2)):9.2f} 1/TeV  '
          f'median |fit-proxy|={np.median(np.abs(df-dp)):8.2f}  frac(|fit-proxy|>14/TeV)={np.mean(np.abs(df-dp)>14.):.3f}  corr={np.corrcoef(dp,df)[0,1]:.3f}')
