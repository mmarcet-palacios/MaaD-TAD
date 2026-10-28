"""extract_fig1b.py -- cut the PRAMEF window out of the H1 Micro-C map for Fig 1B (GRCh38).
Run in MaaD-TAD\\control_loci (where 4DNFI9GMP2J8.mcool lives):
    python extract_fig1b.py
Writes fig1b_microc_H1.npz (about 10 MB) next to the mcool: raw counts, balanced matrix and
balancing weights at 2 kb (matching SEQSIM), 5 kb and 10 kb for GRCh38 chr1:11,800,000-14,400,000 (locus 12.77-13.42 Mb + ~1 Mb flanks).
Reads the file with h5py only (no cooler needed). Duplicate pixels (a cooler 0.8.3 quirk) keep the max, never summed."""
import sys, numpy as np, h5py
MCOOL = sys.argv[1] if len(sys.argv) > 1 else '4DNFI9GMP2J8.mcool'
CHROM, START, END = 'chr1', 11_800_000, 14_400_000
LOCUS = (12_771_346, 13_423_345)

def cut(grp, res):
    names = [c.decode() if isinstance(c, bytes) else c for c in grp['chroms/name'][:]]
    cname = next(n for n in names if n in (CHROM, CHROM[3:])); cid = names.index(cname)
    bc, bs = grp['bins/chrom'][:], grp['bins/start'][:]
    w = grp['bins/weight'][:] if 'weight' in grp['bins'] else np.ones(len(bs))
    sel = np.where((bc == cid) & (bs >= START) & (bs < END))[0]
    lo, hi = sel[0], sel[-1] + 1; n = hi - lo
    off = grp['indexes/bin1_offset']; p0, p1 = int(off[lo]), int(off[hi])
    b1 = grp['pixels/bin1_id'][p0:p1]; b2 = grp['pixels/bin2_id'][p0:p1]; c = grp['pixels/count'][p0:p1]
    k = (b2 >= lo) & (b2 < hi); b1, b2, c = b1[k] - lo, b2[k] - lo, c[k].astype(np.float64)
    raw = np.zeros((n, n)); np.maximum.at(raw, (b1, b2), c); raw = np.maximum(raw, raw.T)
    ww = w[lo:hi].astype(float)
    bal = raw * ww[:, None] * ww[None, :]; bal[np.isnan(ww), :] = np.nan; bal[:, np.isnan(ww)] = np.nan
    return raw.astype(np.float32), bal.astype(np.float32), ww.astype(np.float32), bs[lo:hi].astype(np.int64)

f = h5py.File(MCOOL, 'r')
have = sorted(int(r) for r in f['resolutions'])
out = {'chrom': CHROM, 'window': np.array([START, END]), 'locus': np.array(LOCUS), 'source': MCOOL}
for res in (2000, 5000, 10000):
    if res not in have:
        print(f'{res} bp not in file (have {have}); skipped'); continue
    raw, bal, w, starts = cut(f[f'resolutions/{res}'], res)
    out.update({f'raw_{res}': raw, f'bal_{res}': bal, f'w_{res}': w, f'starts_{res}': starts})
    inl = (starts >= LOCUS[0] - res) & (starts < LOCUS[1])
    cov = np.nansum(raw, 1)
    print(f'{res // 1000} kb: {len(starts)} bins | locus bins masked by balancing: {np.isnan(w[inl]).sum()}/{inl.sum()} '
          f'| median raw coverage: locus {np.median(cov[inl]):.0f} vs flanks {np.median(cov[~inl]):.0f}')
np.savez_compressed('fig1b_microc_H1.npz', **out)
print('wrote fig1b_microc_H1.npz')
