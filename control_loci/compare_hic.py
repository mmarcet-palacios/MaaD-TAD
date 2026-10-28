"""compare_hic.py -- test homology-as-proximity against real Hi-C / Micro-C for one locus.

Reads a .mcool (4DN / cooler format, via h5py only) or a .hic (needs `pip install hic-straw`),
extracts the locus, and compares it with the homology (SEQSIM-like) matrix and the
CSynth-like 3D path built from homology.

  python compare_hic.py --map H1_MicroC.mcool --locus PCDH --res 2000 --fold
  python compare_hic.py --map GM12878.hic   --locus ZNF_chr19 --res 5000

Outputs <locus>_hic_compare.json and <locus>_hic_compare.png next to the inputs.
Tests:
  1 homologous bin pairs vs distance-matched non-homologous pairs (observed/expected contact)
  2 Spearman of Hi-C O/E vs model-predicted O/E (contacts ~ d^-3 from the homology path)
  3 (--fold) fold the Hi-C itself with the same CSynth-like model; compare the two 3D paths
"""
import argparse, json, os, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
LOCI = {  # name: (chrom, start, n_bins at 2 kb) -- GRCh38, same windows as the *_seqsim_like.csv files
    'PRAMEF': ('chr1', 12771321, 325), 'PCDH': ('chr5', 140700000, 450), 'IL13_Th2': ('chr5', 132300000, 350),
    'HBB_OR': ('chr11', 5000000, 350), 'ZNF_chr19': ('chr19', 57000000, 500)}

def read_mcool(path, chrom, start, end, res):
    import h5py
    f = h5py.File(path, 'r')
    if 'resolutions' in f:
        have = sorted(int(r) for r in f['resolutions'])
        if str(res) not in f['resolutions']:
            raise SystemExit(f'resolution {res} not in this file; available: {have}. Rerun with --res <one of these, a multiple of 2000>')
        grp = f[f'resolutions/{res}']
    else:
        grp = f
    names = [c.decode() if isinstance(c, bytes) else c for c in grp['chroms/name'][:]]
    alias = {chrom, chrom.replace('chr', ''), 'chr' + chrom.replace('chr', '')}
    cname = next(n for n in names if n in alias); cid = names.index(cname)
    b_chrom = grp['bins/chrom'][:]; b_start = grp['bins/start'][:]
    w = grp['bins/weight'][:] if 'weight' in grp['bins'] else np.ones(len(b_start))
    sel = np.where((b_chrom == cid) & (b_start >= start - res + 1) & (b_start < end))[0]
    lo, hi = sel[0], sel[-1] + 1; n = hi - lo
    off = grp['indexes/bin1_offset']; p0, p1 = off[lo], off[hi]
    b1 = grp['pixels/bin1_id'][p0:p1]; b2 = grp['pixels/bin2_id'][p0:p1]; c = grp['pixels/count'][p0:p1]
    keep = (b2 >= lo) & (b2 < hi); b1, b2, c = b1[keep] - lo, b2[keep] - lo, c[keep].astype(float)
    # 4DN files made with cooler 0.8.3 can list a pixel more than once; keep the largest count, never sum
    M = np.zeros((n, n)); np.maximum.at(M, (b1, b2), c); M = np.maximum(M, M.T)
    ww = w[lo:hi]; M = M * ww[:, None] * ww[None, :]
    M[np.isnan(ww), :] = np.nan; M[:, np.isnan(ww)] = np.nan
    return M, b_start[lo]

def read_hic(path, chrom, start, end, res):
    import hicstraw
    hic = hicstraw.HiCFile(path); c = chrom.replace('chr', '') if chrom not in [x.name for x in hic.getChromosomes()] else chrom
    mzd = hic.getMatrixZoomData(c, c, 'observed', 'KR', 'BP', res)
    M = mzd.getRecordsAsMatrix(start, end - 1, start, end - 1).astype(float)
    return M, start

def oe(M):
    O = np.full_like(M, np.nan); n = len(M)
    for d in range(n):
        v = np.diagonal(M, d); e = np.nanmean(v) if np.isfinite(v).any() else np.nan
        idx = np.arange(n - d); O[idx, idx + d] = v / e if e else np.nan; O[idx + d, idx] = O[idx, idx + d]
    return O

def coarsen(A, f, fn=np.nanmean):
    n = len(A) // f; A = A[:n * f, :n * f]
    return fn(fn(A.reshape(n, f, n * f), axis=1).reshape(n, n, f), axis=2)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--map', required=True); ap.add_argument('--locus', required=True, choices=LOCI)
    ap.add_argument('--res', type=int, default=2000); ap.add_argument('--dir', default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument('--fold', action='store_true', help='also fold the Hi-C with the CSynth-like model and compare paths')
    a = ap.parse_args()
    chrom, start, n2 = LOCI[a.locus]; end = start + n2 * 2000
    H, h0 = (read_mcool if a.map.endswith(('.mcool', '.cool')) else read_hic)(a.map, chrom, start, end, a.res)
    f = a.res // 2000
    S = np.loadtxt(os.path.join(a.dir, f'{a.locus}_seqsim_like.csv'), delimiter=',', skiprows=1,
                   converters={0: lambda s: 0}, usecols=range(1, n2 + 1))
    X = np.loadtxt(os.path.join(a.dir, f'{a.locus}_csynth_like.xyz'))[:, 1:]
    D = np.linalg.norm(X[:, None] - X[None], axis=-1); P = np.where(D > 0, D, np.nan) ** -3.0
    if f > 1: S = coarsen(S, f, np.nanmax); P = coarsen(P, f)
    n = min(len(H), len(S)); H, S, P = H[:n, :n], S[:n, :n], P[:n, :n]
    Ho, Po = oe(H), oe(P)
    iu = np.triu_indices(n, 3); hv, pv, sv = Ho[iu], Po[iu], S[iu]; ok = np.isfinite(hv) & np.isfinite(pv) & (hv > 0)
    out = dict(locus=a.locus, map=os.path.basename(a.map), res=a.res, bins=n)
    # test 1: homologous vs distance-matched
    hom = ok & (sv > 1.1); rng = np.random.default_rng(0); ratios = []; hv_h = []; hv_c = []   # 1.1 = masked-alignment score 30, the homology threshold used for the loci
    sep = iu[1] - iu[0]
    for k in np.where(hom)[0]:
        cand = np.where(ok & (sv <= 0.6) & (sep == sep[k]))[0]
        if len(cand): c = rng.choice(cand); hv_h.append(hv[k]); hv_c.append(hv[c]); ratios.append(hv[k] / hv[c])
    if ratios:
        out['homologous_pairs'] = len(ratios); out['median_OE_homologous'] = float(np.median(hv_h))
        out['median_OE_matched_controls'] = float(np.median(hv_c)); out['median_ratio'] = float(np.median(ratios))
        out['mannwhitney_p_one_sided'] = float(mannwhitneyu(hv_h, hv_c, alternative='greater').pvalue)
        np.savetxt(os.path.join(a.dir, f'{a.locus}_pairs_OE.csv'), np.c_[hv_h, hv_c], delimiter=',',
                   header='OE_homologous,OE_matched_control', comments='', fmt='%.5f')   # per-pair values for the final figure
    # test 2: model vs Hi-C
    out['spearman_OE_model_vs_hic'] = float(spearmanr(np.log(pv[ok]), np.log(hv[ok])).correlation)
    # test 3: fold the Hi-C and compare paths
    if a.fold and f == 1:
        from csynth_like import fold
        Hs = np.where(np.isfinite(H), H, 0); Hs = 100 * Hs / np.nanpercentile(Hs[iu], 99.9); np.fill_diagonal(Hs, 100)
        Xh, _ = fold(np.clip(Hs, 0, 100), kc=0.25, rep=0.3, seeds=4)
        Dh = np.linalg.norm(Xh[:, None] - Xh[None], axis=-1)
        out['spearman_paths_homology_vs_hic'] = float(spearmanr(D[:n, :n][iu], Dh[iu]).correlation)
        np.savetxt(os.path.join(a.dir, f'{a.locus}_hic_fold.xyz'), np.c_[start + np.arange(n) * 2000, Xh], fmt=['%d', '%.3f', '%.3f', '%.3f'], delimiter='\t')
    json.dump(out, open(os.path.join(a.dir, f'{a.locus}_hic_compare.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        for x, M, t in zip(ax, [np.log2(Ho), np.log(S), np.log2(Po)], ['Hi-C  log2(O/E)', 'homology (SEQSIM-like, log)', 'homology 3D path  log2(O/E)']):
            x.imshow(M, cmap='RdBu_r' if 'O/E' in t else 'Oranges', interpolation='nearest'); x.set_title(t); x.set_xticks([]); x.set_yticks([])
        fig.suptitle(f"{a.locus} · {out['map']} · {a.res // 1000} kb"); plt.savefig(os.path.join(a.dir, f'{a.locus}_hic_compare.png'), bbox_inches='tight')
    except Exception as e: print('plot skipped:', e)
