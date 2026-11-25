"""analyse_ensemble.py -- what is shared, and what varies, across the MaaD-TAD path ensemble.

Reads ensemble/m*/ (from ensemble_member.py) and the published model (data/, frames_tethered.npz), writes
ensemble/summary.json, ensemble/members.csv, ensemble/blocks.csv and ensemble/maps.npz.

Per member: domain size (radius of gyration and maximum extent of the nucleosome centres); homologous-pair contact
enrichment (median d^-3 ratio, homologous vs same-separation control; as in Fig 8C) at each Folder stage; pairing of
every homology block at nucleosome level (spacing, CV, direction, register; as pairing_scorecard.py); tether
distances; homologous gene-pair distances; nucleosome collisions.
Across members: path-to-path similarity (Spearman correlation of the 2 kb bead distance matrices) and the
ensemble contact map, compared with H1 Micro-C (Spearman of observed/expected, single members vs ensemble mean).
"""
import glob, json, os, shutil, sys
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.abspath(__file__))
EN = os.path.join(ROOT, 'ensemble')
N = 325
LO = 12771273          # ROI_START of stage1_frames.py; bead k = section k + 163
GAPS = {int(l) for l in open(os.path.join(ROOT, 'data', 'gap_sections.txt')) if l.strip()}
R_ALL = json.load(open(os.path.join(ROOT, 'data', 'homology_restraints.json')))
TETH = {t['id']: (0.5 * (t['cohesin'][0] + t['cohesin'][1]), t['mirror']) for t in R_ALL['tethers']}

# homology map and control pairs (as pramef_calibration.py)
S = np.loadtxt(os.path.join(ROOT, 'control_loci', 'PRAMEF_seqsim_like.csv'), delimiter=',', skiprows=1,
               converters={0: lambda s: 0}, usecols=range(1, N + 1))[:N, :N]
iu = np.triu_indices(N, 3); sep = iu[1] - iu[0]; sv = S[iu]; hom = sv > 1.1
rng = np.random.default_rng(0); hk = np.where(hom)[0]
ctrl = np.array([rng.choice(np.where((sv <= 0.6) & (sep == sep[i]))[0]) for i in hk])
real = np.array([(k + 163) not in GAPS for k in range(N)])   # gap_sections.txt uses section numbers from 163


def enr(X):
    D = np.linalg.norm(X[:, None] - X[None], axis=-1); P = D[iu] ** -3.0
    return float(np.median(P[hk] / P[ctrl]))


def xyz(p):
    return np.loadtxt(p)[:N, 1:]


def scale_like_stage1(X, domain_nm):
    r = X[real]; dmax = np.linalg.norm(r[:, None] - r[None], axis=-1).max()
    return X * (domain_nm * 10.0 / dmax)


# homologous gene pairs: a gene on arm A and the gene whose centre lies at the homologous position on arm C
GENES = [l.split('\t') for l in open(os.path.join(ROOT, 'figures', 'data', 'genes.tsv')) if l.strip()]
GENES = [(g[0], 0.5 * (int(g[2]) + int(g[3]))) for g in GENES]


def partner(b, g):
    A0, A1 = b['armA']; C0, C1 = b['armC']; t = (g - A0) / (A1 - A0)
    return C0 + t * (C1 - C0) if b['type'] == 'direct' else C1 - t * (C1 - C0)


GPAIRS = []
for b in R_ALL['blocks']:
    for name, g in GENES:
        if b['armA'][0] <= g <= b['armA'][1]:
            p = partner(b, g)
            cand = [(abs(gc - p), n2, gc) for n2, gc in GENES if b['armC'][0] <= gc <= b['armC'][1]]
            if cand and min(cand)[0] < 5000:
                GPAIRS.append((b['id'], name, min(cand)[1], g, min(cand)[2]))


def frame_metrics(F, tethers_on, pinch_kb=6.0):
    O = F['origins']; dy = F['dyad_pos'].astype(float)
    pos = lambda g: np.column_stack([np.interp(g, dy, O[:, k]) for k in range(3)])
    c = O - O.mean(0)
    out = dict(n_nuc=len(O), rg_nm=float(np.sqrt((c ** 2).sum(1).mean()) / 10),
               extent_nm=float(np.linalg.norm(O[:, None][::7] - O[None][:, ::7], axis=-1).max() / 10))
    T = [(tid, *TETH[tid]) for tid in TETH]   # same samples in every member: skip all three tether zones
    blocks = {}
    for b in R_ALL['blocks']:
        A0, A1 = b['armA']; C0, C1 = b['armC']
        n = int((A1 - A0) / 1000) + 1; t = np.linspace(0, 1, n)
        gA = A0 + t * (A1 - A0); gC = C0 + t * (C1 - C0) if b['type'] == 'direct' else C1 - t * (C1 - C0)
        keep = np.ones(n, bool)
        for _, gc, gm in T:
            for g in (gA, gC):
                keep &= np.minimum(np.abs(g - gc), np.abs(g - gm)) > pinch_kb * 1000
        if keep.sum() < 3:
            continue
        PA, PC = pos(gA), pos(gC); d = np.linalg.norm(PA - PC, axis=1)[keep]
        h = 1000.0; tA = pos(gA + h) - pos(gA - h); tC = pos(gC + h) - pos(gC - h)
        tA /= np.linalg.norm(tA, axis=1)[:, None]; tC /= np.linalg.norm(tC, axis=1)[:, None]
        sg = 1 if b['type'] == 'direct' else -1
        ang = np.degrees(np.arccos(np.clip(sg * np.sum(tA * tC, 1), -1, 1)))[keep]
        gCd = np.arange(C0, C1 + 1, 50.0); near = gCd[cKDTree(pos(gCd)).query(PA[keep])[1]]
        blocks[b['id']] = dict(spacing_nm=float(np.median(d) / 10), cv=float(d.std() / d.mean()),
                               dir_deg=float(np.median(ang)), register_bp=float(np.median(np.abs(near - gC[keep]))))
    out['blocks'] = blocks
    out['tether_nm'] = {tid: float(np.linalg.norm(pos([gc]) - pos([gm])) / 10) for tid, (gc, gm) in TETH.items()}
    out['genepair_nm'] = {f'{a}-{b2}': float(np.linalg.norm(pos([ga]) - pos([gb])) / 10) for _, a, b2, ga, gb in GPAIRS}
    p = cKDTree(O).query_pairs(110.0, output_type='ndarray'); p = p[np.abs(p[:, 0] - p[:, 1]) > 2]
    out['collisions'] = int(len(p))
    # 2 kb bead positions from the nucleosome centres (for maps and path similarity)
    out['_beads'] = pos(LO + 1000 + 2000 * np.arange(N))
    return out


def stage1_log(d):
    import re
    p = os.path.join(d, 'log.txt')
    if not os.path.exists(p):
        return {}
    t = open(p).read().split('> stage1_frames.py')[-1].split('\n> ')[0]
    g = lambda pat: re.search(pat, t)
    out = {}
    m = g(r'over 1.0 \(unreachable\): (\d+)/(\d+)'); out['unreachable'] = int(m.group(1)) if m else None
    m = g(r'aim\s+: median ([0-9.]+) deg, over 28: (\d+)/(\d+)'); out['aim_over28_frac'] = int(m.group(2)) / int(m.group(3)) if m else None
    m = g(r'steering correction : median ([0-9.]+) deg/nuc, 95th ([0-9.]+)'); out['steer_median'] = float(m.group(1)) if m else None
    m = g(r'extension\s+: median ([0-9.]+)'); out['ext_median'] = float(m.group(1)) if m else None
    return out


def load_member(d, info=None):
    m = dict(name=os.path.basename(d))
    if info is None:
        info = json.load(open(os.path.join(d, 'member.json')))
    m['tethers'] = info['tethers']; m['seed'] = info.get('seed')
    m['median_nrl'] = info.get('median_nrl')
    fr = os.path.join(d, 'frames.npz') if os.path.exists(os.path.join(d, 'frames.npz')) and info.get('times', {}).get('relax') else os.path.join(d, 'frames_stage1.npz')
    m['frames_file'] = os.path.basename(fr)
    m.update(frame_metrics(np.load(fr), m['tethers']))
    m['stage1'] = stage1_log(d)
    from register_path import register, pairing
    m['pairing_path'] = {'refolded': pairing(os.path.join(d, 'guide_refolded.xyz'), info.get('domain_refold', 963.5)),
                         'tethered': pairing(os.path.join(d, 'guide_tethered.xyz'), info.get('domain_tether', 1011.5))}
    m['register_path'] = {st: register(os.path.join(d, f)) for st, f in (('fold', 'fold.xyz'), ('refolded', 'guide_refolded.xyz'), ('tethered', 'guide_tethered.xyz'))}
    m['enr'] = dict(fold=enr(xyz(os.path.join(d, 'fold.xyz'))), refolded=enr(xyz(os.path.join(d, 'guide_refolded.xyz'))),
                    tethered=enr(xyz(os.path.join(d, 'guide_tethered.xyz'))), nucleosomes=enr(m['_beads']))
    return m


def oe(M):
    M = M.copy(); n = len(M)
    for k in range(n):
        idx = np.arange(n - k); v = M[idx, idx + k]; mu = np.nanmean(v)
        M[idx, idx + k] = v / mu; M[idx + k, idx] = v / mu
    return M


def main():
    global EN
    if len(sys.argv) > 1:
        EN = os.path.join(ROOT, sys.argv[1])
    dirs = sorted(d for d in glob.glob(os.path.join(EN, 'm*')) if os.path.exists(os.path.join(d, 'member.json')))
    mem = [load_member(d) for d in dirs]
    pubd = os.path.join(EN, '_published'); os.makedirs(pubd, exist_ok=True)
    for src, dst in (('gap-filled.xyz', 'fold.xyz'), ('guide_refolded.xyz', 'guide_refolded.xyz'), ('guide_tethered.xyz', 'guide_tethered.xyz')):
        if not os.path.exists(os.path.join(pubd, dst)):
            shutil.copy(os.path.join(ROOT, 'data', src), os.path.join(pubd, dst))   # copies, not symlinks (Windows)
    if not os.path.exists(os.path.join(pubd, 'frames_stage1.npz')):
        shutil.copy(os.path.join(ROOT, 'frames_tethered.npz'), os.path.join(pubd, 'frames_stage1.npz'))
    pub = load_member(pubd, dict(tethers=list(TETH), seed=None, median_nrl=190.0, domain_refold=963.5, domain_tether=1011.5))
    pub['name'] = 'published'
    allm = mem + [pub]
    # path-to-path similarity on 2 kb beads from nucleosome centres (gap sections excluded)
    B = [m.pop('_beads')[real] for m in allm]
    Ds = [np.linalg.norm(b[:, None] - b[None], axis=-1) for b in B]
    ju = np.triu_indices(real.sum(), 3)
    k = len(allm); rho = np.eye(k)
    for i in range(k):
        for j in range(i + 1, k):
            rho[i, j] = rho[j, i] = spearmanr(Ds[i][ju], Ds[j][ju]).correlation
    # contact maps: d^-3 per member (beads), ensemble mean over members (not the published one)
    P = np.array([np.where(D > 0, D, np.nan) ** -3.0 for D in Ds])
    Pm = np.nanmean(P[:-1], axis=0)
    # H1 Micro-C at 2 kb, same beads
    z = np.load(os.path.join(ROOT, 'control_loci', 'fig1b_microc_H1.npz'))
    st = z['starts_2000']; H = z['bal_2000'].astype(float)
    sel = np.array([int(np.argmin(np.abs(st + 1000 - (LO + 1000 + 2000 * k2)))) for k2 in range(N)])[real]
    Hm = H[np.ix_(sel, sel)]
    Ho = oe(Hm); ok = np.isfinite(Ho[ju]) & (Ho[ju] > 0)
    def vs_hic(M):
        Mo = oe(M)[ju]; g = ok & np.isfinite(Mo)
        return float(spearmanr(np.log(Mo[g]), np.log(Ho[ju][g])).correlation)
    hic_single = [vs_hic(Pk) for Pk in P]
    hic_ens = vs_hic(Pm)
    # ensemble-mean homologous enrichment: on the mean map
    Pfull = np.full((N, N), np.nan); Pfull[np.ix_(real, real)] = Pm
    ens_enr = float(np.nanmedian(Pfull[iu][hk] / Pfull[iu][ctrl]))
    # running mean: how fast does the ensemble map converge?
    conv = []
    for n_ in (1, 2, 5, 10, 20, 30, 40, 50):
        if n_ <= len(mem):
            conv.append((n_, vs_hic(np.nanmean(P[:n_], axis=0))))
    # summaries
    def stats(v):
        v = np.array([x for x in v if x is not None and np.isfinite(x)], float)
        return dict(n=int(len(v)), median=float(np.median(v)), p5=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)),
                    min=float(v.min()), max=float(v.max())) if len(v) else None
    S_ = dict(n_members=len(mem), gene_pairs=[g[:3] for g in GPAIRS])
    for key in ('rg_nm', 'extent_nm', 'n_nuc', 'collisions', 'median_nrl'):
        S_[key] = stats([m[key] for m in mem]); S_[key]['published'] = pub[key]
    for st_ in ('fold', 'refolded', 'tethered', 'nucleosomes'):
        S_['enr_' + st_] = stats([m['enr'][st_] for m in mem]); S_['enr_' + st_]['published'] = pub['enr'][st_]
    S_['enr_ensemble_mean_map'] = ens_enr
    r_mm = rho[:-1, :-1][np.triu_indices(len(mem), 1)]; r_mp = rho[-1, :-1]
    S_['path_similarity'] = dict(member_member=stats(r_mm), member_published=stats(r_mp))
    S_['hic_spearman'] = dict(single=stats(hic_single[:-1]), published=hic_single[-1], ensemble_mean=hic_ens, convergence=conv)
    for st in ('fold', 'refolded', 'tethered'):
        S_['register_' + st] = stats([float(np.median(list(m['register_path'][st].values()))) for m in mem])
        S_['register_' + st]['published'] = float(np.median(list(pub['register_path'][st].values())))
    for st in ('refolded', 'tethered'):
        for q in ('slide_nm', 'dir_deg', 'spacing_nm', 'cv'):
            S_[f'pair_{st}_{q}'] = stats([float(np.median([v[q] for v in m['pairing_path'][st].values()])) for m in mem])
            S_[f'pair_{st}_{q}']['published'] = float(np.median([v[q] for v in pub['pairing_path'][st].values()]))
    S_['pair_block_tethered'] = {b: {q: stats([m['pairing_path']['tethered'][b][q] for m in mem]) for q in ('slide_nm', 'dir_deg', 'spacing_nm')}
                                 for b in mem[0]['pairing_path']['tethered']}
    S_['register_block_tethered'] = {b: stats([m['register_path']['tethered'][b] for m in mem]) for b in mem[0]['register_path']['tethered']}
    for q in ('unreachable', 'aim_over28_frac', 'steer_median', 'ext_median'):
        S_['stage1_' + q] = stats([m['stage1'].get(q) for m in mem])
    S_['blocks'] = {}
    for b in R_ALL['blocks']:
        rows = [m['blocks'][b['id']] for m in mem if b['id'] in m['blocks']]
        if not rows:
            continue
        S_['blocks'][b['id']] = {q: stats([r[q] for r in rows]) for q in ('spacing_nm', 'cv', 'dir_deg', 'register_bp')}
        S_['blocks'][b['id']]['published'] = pub['blocks'].get(b['id'])
        S_['blocks'][b['id']]['type'] = b['type']
    S_['tethers'] = {}
    for tid in TETH:
        on = [m['tether_nm'][tid] for m in mem if tid in m['tethers']]
        off = [m['tether_nm'][tid] for m in mem if tid not in m['tethers']]
        S_['tethers'][tid] = dict(on=stats(on), off=stats(off), published=pub['tether_nm'][tid])
    S_['genepairs'] = {}
    for gp in mem[0]['genepair_nm']:
        S_['genepairs'][gp] = stats([m['genepair_nm'][gp] for m in mem]); S_['genepairs'][gp]['published'] = pub['genepair_nm'][gp]
    json.dump(S_, open(os.path.join(EN, 'summary.json'), 'w'), indent=1)
    # tables
    with open(os.path.join(EN, 'members.csv'), 'w') as f:
        f.write('member,frames,tethers,n_nuc,median_nrl,rg_nm,extent_nm,collisions,enr_fold,enr_refolded,enr_tethered,enr_nucleosomes,'
                'rho_to_published,hic_spearman,' + ','.join(f'{t}_nm' for t in TETH) + '\n')
        for i, m in enumerate(allm):
            f.write(','.join(str(x) for x in [m['name'], m['frames_file'], '+'.join(m['tethers']) or 'none', m['n_nuc'], m['median_nrl'],
                    round(m['rg_nm'], 1), round(m['extent_nm'], 1), m['collisions']] + [round(m['enr'][s], 2) for s in ('fold', 'refolded', 'tethered', 'nucleosomes')]
                    + [round(rho[i, -1], 3), round(hic_single[i], 3)] + [round(m['tether_nm'][t], 1) for t in TETH]) + '\n')
    with open(os.path.join(EN, 'blocks.csv'), 'w') as f:
        f.write('member,block,spacing_nm,cv,dir_deg,register_bp\n')
        for m in allm:
            for bid, r in m['blocks'].items():
                f.write(f"{m['name']},{bid},{r['spacing_nm']:.1f},{r['cv']:.3f},{r['dir_deg']:.1f},{r['register_bp']:.0f}\n")
    with open(os.path.join(EN, 'pairing_blocks.csv'), 'w') as f:
        f.write('member,stage,block,spacing_nm,cv,dir_deg,slide_nm\n')
        for m in allm:
            for st in ('refolded', 'tethered'):
                for bid, r in m['pairing_path'][st].items():
                    f.write(f"{m['name']},{st},{bid},{r['spacing_nm']:.1f},{r['cv']:.3f},{r['dir_deg']:.1f},{r['slide_nm']:.2f}\n")
    np.savez_compressed(os.path.join(EN, 'maps.npz'), rho=rho, names=[m['name'] for m in allm], P_mean=Pm, P_pub=P[-1],
                        P_members=P[:-1].astype(np.float32), hic=Hm, real=real, hic_single=hic_single)
    print(json.dumps({k: v for k, v in S_.items() if k not in ('blocks', 'genepairs', 'gene_pairs')}, indent=1))


if __name__ == '__main__':
    main()
