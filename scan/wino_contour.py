"""Compare wino_validation_{chgchg,chgn1}.csv with (i) the public recast efficiencies and (ii) the official ATLAS contour.
Exclusion: sigma_NLO+NLL(chgchg + chgn1)(m) x eff_total > 0.037 fb (ATLAS-SUSY-2018-19 model-independent sigma_vis limit).  usage: wino_contour.py <ATLAS_data> <data dir>"""
import sys, numpy as np, pandas as pd
atl, dd = sys.argv[1:3]
d = pd.concat([pd.read_csv(f'{dd}/wino_validation_{p}.csv') for p in ('chgchg', 'chgn1')])
d['ratio'] = d.eff / d.eff_public_recast
print('efficiency / public recast efficiency (Pythia-SUSY winos, generator-level chain):')
print(d.pivot_table(index=['proc', 'M', 'tau_ns'], columns='smear', values='ratio').round(2).to_string())
xs = {p: pd.read_csv(f'{dd}/xsec_13TeV_pp2{p}_fb.csv').set_index('mass').xsec for p in ('chgchg', 'chgn1')}
pts = np.array([tuple(map(float, l.split(';'))) for l in open(f'{atl}/official_exclusion.txt') if ';' in l])
low, up = pts[:49], pts[48:]            # lower-tau branch (m rising with tau) and upper-tau branch (m falling with tau)
def official_mass(tau):
    # largest mass excluded at this lifetime: where tau crosses the lower branch (tau_min(m) = tau) or the upper branch
    m_low = np.interp(tau, low[:, 1], low[:, 0]) if low[:, 1].min() < tau < low[:, 1].max() else np.nan
    o = np.argsort(up[:, 1]); m_up = np.interp(tau, up[o, 1], up[o, 0]) if up[:, 1].min() < tau < up[:, 1].max() else np.nan
    return m_low, m_up
print('\nmass limit (GeV) from sigma x eff > 0.037 fb, interpolated in M (log):')
for tau in sorted(d.tau_ns.unique()):
    out = {}
    for conv in ('first', 'intended'):
        Ms, v = [], []
        for M in sorted(d.M.unique()):
            g = d[(d.M == M) & (d.tau_ns == tau) & (d.smear == conv)].set_index('proc')
            if len(g) < 2: continue
            s = {p: float(np.interp(M, xs[p].index, xs[p].values)) for p in xs}
            Ms.append(M); v.append((s['chgchg'] * g.loc['chgchg', 'eff'] + s['chgn1'] * g.loc['chgn1', 'eff']))
        Ms, v = np.array(Ms, float), np.array(v)
        ok = np.where(v > 0.037)[0]
        if len(ok) == 0: out[conv] = '<%d' % Ms.min(); continue
        k = ok[-1]
        if k + 1 < len(Ms): out[conv] = '%.0f' % np.exp(np.interp(np.log(0.037), [np.log(v[k + 1]), np.log(v[k])], [np.log(Ms[k + 1]), np.log(Ms[k])]))
        else: out[conv] = '>%d' % Ms[-1]
    ml, mu = official_mass(tau)
    print(f'tau={tau:4.1f} ns: first {out["first"]:>6}  intended {out["intended"]:>6}   official contour: lower-branch m={ml:.0f}  upper-branch m={mu:.0f}')
