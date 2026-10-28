"""sweep_calibrate.py -- calibrate the homology weight against Micro-C (Fig 9C).

For each locus and homology weight k, fold the homology map (csynth_like.fold) from several random starts
and measure (1) the contact enrichment of homologous vs same-distance control pairs in the model
(contacts ~ d^-3) and (2) agreement of the path with the Micro-C fold (<locus>_hic_fold.xyz, from
compare_hic.py --fold). The calibrated k is where (1) matches the Micro-C enrichment.

  python sweep_calibrate.py            # uses all CPU cores; ~10 min on 24 cores, ~1 h on 2
Writes sweep_calibrate.json and sweep_calibrate.png.
"""
import json, os, numpy as np
from multiprocessing import Pool
from scipy.stats import spearmanr
from csynth_like import fold
HERE = os.path.dirname(os.path.abspath(__file__))
LOCI = {'PCDH': (450, None), 'ZNF_chr19': (500, None)}
KS = [0, 0.002, 0.005, 0.01, 0.02, 0.05]
MODES = {'homologous pairs only': 1.1, 'whole map': None}   # threshold on the SEQSIM-like scale
SEEDS = 5

def load(loc):
    n = LOCI[loc][0]
    S = np.loadtxt(os.path.join(HERE, f'{loc}_seqsim_like.csv'), delimiter=',', skiprows=1,
                   converters={0: lambda s: 0}, usecols=range(1, n + 1))
    H = np.loadtxt(os.path.join(HERE, f'{loc}_hic_fold.xyz'))[:, 1:]
    obs = json.load(open(os.path.join(HERE, f'{loc}_hic_compare.json')))['median_ratio']
    return S, H, obs

def job(args):
    loc, mode, k, seed = args
    S, H, obs = load(loc); n = len(H); S = S[:n, :n]
    iu = np.triu_indices(n, 3); sep = iu[1] - iu[0]; sv = S[iu]
    hom = sv > 1.1; rng = np.random.default_rng(0); hk = np.where(hom)[0]
    ctrl = np.array([rng.choice(np.where((sv <= 0.6) & (sep == sep[i]))[0]) for i in hk])
    M = S if k > 0 else np.full((n, n), 0.4)
    X, _ = fold(M, kc=(k if k > 0 else 0.25), rep=0.3, seeds=1, rng0=100 + seed, thresh=MODES[mode])
    D = np.linalg.norm(X[:, None] - X[None], axis=-1); P = D[iu] ** -3.0
    Dh = np.linalg.norm(H[:, None] - H[None], axis=-1)
    # agreement beyond genomic distance: z-score log distances within each separation, then correlate
    def zsep(v):
        v = np.log(v); z = np.empty_like(v)
        for d in np.unique(sep):
            m = sep == d; sd = v[m].std(); z[m] = (v[m] - v[m].mean()) / (sd if sd > 0 else 1)
        return z
    beyond = float(np.corrcoef(zsep(D[iu]), zsep(Dh[iu]))[0, 1])
    return loc, mode, k, seed, float(np.median(P[hk] / P[ctrl])), float(spearmanr(D[iu], Dh[iu]).correlation), obs, beyond

if __name__ == '__main__':
    tasks = [(l, m, k, s) for l in LOCI for m in MODES for k in KS for s in range(SEEDS)]
    with Pool(max(1, (os.cpu_count() or 2) - 2) if (os.cpu_count() or 2) > 4 else (os.cpu_count() or 1)) as p:
        rows = []
        for i, r in enumerate(p.imap_unordered(job, tasks)):
            rows.append(r); print(f'[{i + 1}/{len(tasks)}] {r[0]} {r[1]} k={r[2]} seed {r[3]}: enrichment {r[4]:.2f}', flush=True)
    out = {}
    for loc in LOCI:
        obs = [r[6] for r in rows if r[0] == loc][0]
        out[loc] = {'observed_enrichment': obs}
        for m in MODES:
            ws = []
            for k in KS:
                e = [r[4] for r in rows if r[0] == loc and r[1] == m and r[2] == k]
                a = [r[5] for r in rows if r[0] == loc and r[1] == m and r[2] == k]
                b = [r[7] for r in rows if r[0] == loc and r[1] == m and r[2] == k]
                ws.append(dict(k=k, enrichment=float(np.mean(e)), enrichment_sd=float(np.std(e)),
                               agreement=float(np.mean(a)), agreement_sd=float(np.std(a)),
                               agreement_beyond_distance=float(np.mean(b)), agreement_beyond_distance_sd=float(np.std(b))))
                print(f'{loc:10s} {m:22s} k={k:<6} enrichment {np.mean(e):5.2f} ± {np.std(e):.2f} (Micro-C {obs:.2f})   '
                      f'path vs Micro-C fold {np.mean(a):.3f} ± {np.std(a):.3f}   beyond distance {np.mean(b):.3f} ± {np.std(b):.3f}')
            ks = np.array(KS); en = np.maximum.accumulate([w['enrichment'] for w in ws])   # enforce monotone for interpolation
            cal = float(np.interp(obs, en, ks)) if en[0] <= obs <= en[-1] else None
            out[loc][m] = {'weights': ws, 'calibrated_k': cal}
            print(f'{loc} [{m}]: calibrated homology weight k = {cal}')
    json.dump(out, open(os.path.join(HERE, 'sweep_calibrate.json'), 'w'), indent=1)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for loc, c in zip(LOCI, ['#d9731f', '#3f6fa8']):
            for m, ls in zip(MODES, ['-', ':']):
                w = out[loc][m]['weights']; x = [max(v['k'], 0.001) for v in w]
                ax[0].errorbar(x, [v['enrichment'] for v in w], [v['enrichment_sd'] for v in w], color=c, ls=ls, marker='o', label=f'{loc}, {m}')
                ax[1].errorbar(x, [v['agreement_beyond_distance'] for v in w], [v['agreement_beyond_distance_sd'] for v in w], color=c, ls=ls, marker='o', label=f'{loc}, {m}')
            ax[0].axhline(out[loc]['observed_enrichment'], color=c, ls='--', lw=1)
        ax[0].set_yscale('log')
        for a, t in zip(ax, ['model contact enrichment (dashed: Micro-C)', 'path agreement with Micro-C fold, beyond genomic distance']):
            a.set_xscale('log'); a.set_xlabel('homology weight k (0 plotted at 0.001)'); a.set_title(t, fontsize=10); a.legend(fontsize=7)
        plt.tight_layout(); plt.savefig(os.path.join(HERE, 'sweep_calibrate.png'), dpi=150)
    except Exception as e: print('plot skipped:', e)
