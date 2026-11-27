"""Validation of the map-based DT recast chain (generator-level jets/MET + trigger/tracklet maps + pT smearing) on pure winos, against
(i) the efficiencies of the public recast (Ramos & Lessa, validation/C1C1_effs.json, C1N1_effs.json: MG5+Pythia+Delphes, EWK SR, overall) and
(ii) the official ATLAS exclusion contour (ATLAS_data/official_exclusion.txt).
stage 'gen':  wino_validation.py gen <mass> <proc chgchg|chgn1> <nevents> <out.npz>   (Pythia8 SUSY processes with ISR/FSR, hadron level, stable charginos)
stage 'eff':  wino_validation.py eff <ATLAS_data dir> <validation dir> <kfactor csv> <out.csv> <npz dir>   (needs the gen files wino_<proc>_<mass>.npz)
Differences to the public recast: Pythia SUSY processes (no matched ME jets), generator-level jets/MET with the parametrised smearing,
straight-line decay positions."""
import sys, os, json, numpy as np
stage = sys.argv[1]
if stage == 'gen':
    import pythia8, fastjet
    M = float(sys.argv[2]); proc = sys.argv[3]; N = int(sys.argv[4]); out = sys.argv[5]; H = 5.0e4
    src = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'wino_pythia_check.py')).read()
    slha = src[src.index('slha = f"""') + len('slha = f"""'):src.index('"""\nopen(')].replace('{H}', str(H)).replace('{M}', str(M)).replace('{H+10}', str(H + 10)).replace('{H+20}', str(H + 20)).replace('{H*H}', str(H * H))
    open(f'/tmp/wino_{proc}_{int(M)}.spc', 'w').write(slha)
    p = pythia8.Pythia("", False)
    for s in ["Beams:eCM = 13000", "PDF:pSet = 20", "SUSY:all = off", "SUSY:qqbar2chi+chi- = " + ('on' if proc == 'chgchg' else 'off'), "SUSY:qqbar2chi+-chi0 = " + ('on' if proc == 'chgn1' else 'off'),
              f"SLHA:file = /tmp/wino_{proc}_{int(M)}.spc", "SLHA:minMassSM = 100.", "1000024:mayDecay = off", "Next:numberCount = 0", "Random:setSeed = on", f"Random:seed = {int(M) + (7 if proc == 'chgn1' else 3)}",
              "ParticleDecays:limitTau0 = on", "ParticleDecays:tau0Max = 10"]: p.readString(s)
    p.init()
    jd = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    J, MX, MY, CH = [], [], [], []
    for _ in range(N):
        if not p.next(): continue
        vis, mx, my, ch = [], 0., 0., []
        for i in range(p.event.size()):
            pt = p.event[i]
            if not pt.isFinal(): continue
            ida = abs(pt.id())
            if ida in (12, 14, 16, 1000022, 1000024):
                if ida == 1000024: ch.append((pt.px(), pt.py(), pt.pz(), pt.e(), 1. if pt.id() > 0 else -1.))
                continue
            if abs(pt.eta()) < 4.9: mx += pt.px(); my += pt.py()
            if abs(pt.eta()) < 4.7 and pt.pT() > 0.: vis.append(fastjet.PseudoJet(pt.px(), pt.py(), pt.pz(), pt.e()))
        cs = fastjet.ClusterSequence(vis, jd); jets = sorted(cs.inclusive_jets(20.), key=lambda j: -j.pt())[:6]
        a = np.zeros((6, 3))
        for k, j in enumerate(jets): a[k] = (j.pt(), j.eta(), j.phi())
        J.append(a); MX.append(-mx); MY.append(-my); CH.append(np.array(ch))
    np.savez(out, jets=np.array(J), metx=np.array(MX), mety=np.array(MY), ch=np.array(CH, dtype=object), allow_pickle=True)
    print('wrote', out, len(J))
else:
    atl, val, kf, out, ndir = sys.argv[2:7]
    src = open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'atlas_dt_recast.py')).read()
    head = src[:src.index('params = {}')]
    sys.argv = [sys.argv[0], ndir, atl, kf, out, 'first']
    import re, csv, glob, uproot
    exec(head)
    G = 299.792458   # mm/ns
    ref = {}
    for proc, fn in (('chgchg', 'C1C1_effs.json'), ('chgn1', 'C1N1_effs.json')):
        for x in json.load(open(f'{val}/{fn}')):
            for e in x['Efficiencies']: ref[(proc, int(x['mLLP']), round(e['tau_ns'], 4))] = (e['EWK SR'], e['EWK SR Error'])
    taus = [0.05, 0.1, 0.4, 0.8, 2.0, 7.0, 20.0]; rows = []
    rng3 = np.random.default_rng(5)
    for proc in ('chgchg', 'chgn1'):
        for M in (100, 200, 300, 400, 500, 600, 700, 800, 900):
            fn = f'{ndir}/wino_{proc}_{M}.npz'
            if not os.path.exists(fn): continue
            d = np.load(fn, allow_pickle=True); n = len(d['metx']); jm = dict(jets=d['jets'].copy(), metx=d['metx'].copy(), mety=d['mety'].copy())
            # per event objects and trigger/kinematics: once per event (smearing random but shared by all tau)
            obj = []
            for k in range(n):
                pt, eta, ph, met, mph = event_objects(jm, k); trig, kin = presel(pt, eta, ph, met, mph); obj.append((eta, ph, trig if kin else 0.))
            for tau in taus:
                for conv in ('first', 'intended'):
                    SMEAR = conv; tot = 0.; tot2 = 0.
                    for k in range(n):
                        eta, ph, tr = obj[k]
                        if tr <= 0.: continue
                        pr = 1.
                        for c in d['ch'][k]:
                            px, py, pz, e, sq = c; ptq = np.hypot(px, py); p = np.sqrt(px * px + py * py + pz * pz); bg = p / np.sqrt(max(e * e - p * p, 1e-9))
                            l = bg * G * tau * rng3.exponential(); R = l * ptq / p
                            q = dict(pt=ptq, pt0=ptq, eta=np.arcsinh(pz / max(ptq, 1e-9)), phi=np.arctan2(py, px), R=R, dphis=0., L=0.1, sq=sq, u=rng3.random(), chi2ok=1.)
                            pr *= 1. - track_eff(q, eta, ph, False)
                        v = tr * (1. - pr); tot += v; tot2 += v * v
                    eff = tot / n; err = np.sqrt(max(tot2 / n - eff ** 2, 0.) / n)
                    r = ref.get((proc, M, round(tau, 4)), (np.nan, np.nan))
                    rows.append(dict(proc=proc, M=M, tau_ns=tau, smear=conv, eff=eff, err=err, eff_public_recast=r[0], err_public_recast=r[1]))
                    print(proc, M, tau, conv, f'eff={eff:.3e}+-{err:.1e}  public recast {r[0]:.3e}+-{r[1]:.1e}', flush=True)
    with open(out, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
