"""extract_stripe_maps.py -- 5 kb Micro-C windows around the 145 kb D3/D4 line, H1 and HFFc6 (GRCh38).
Run in MaaD-TAD\\control_loci (where both mcool files live):
    python extract_stripe_maps.py
Writes stripe_maps.npz (well under 1 MB): raw counts, balanced matrix, weights and bin starts
for chr1:13,100,000-13,500,000 at 5 kb, one set per file (keys prefixed H1_ and HFFc6_).
h5py only; duplicate pixels keep the max (same rule as extract_fig1b.py)."""
import numpy as np, h5py

FILES = {'H1': '4DNFI9GMP2J8.mcool', 'HFFc6': '4DNFI643OYP9.mcool'}
CHROM, START, END, RES = 'chr1', 13_100_000, 13_500_000, 5000


def cut(grp):
    names = [c.decode() if isinstance(c, bytes) else c for c in grp['chroms/name'][:]]
    cid = names.index(next(n for n in names if n in (CHROM, CHROM[3:])))
    bc, bs = grp['bins/chrom'][:], grp['bins/start'][:]
    w = grp['bins/weight'][:] if 'weight' in grp['bins'] else np.ones(len(bs))
    sel = np.where((bc == cid) & (bs >= START) & (bs < END))[0]
    lo, hi = sel[0], sel[-1] + 1
    off = grp['indexes/bin1_offset']; p0, p1 = int(off[lo]), int(off[hi])
    b1 = grp['pixels/bin1_id'][p0:p1]; b2 = grp['pixels/bin2_id'][p0:p1]; c = grp['pixels/count'][p0:p1]
    k = (b2 >= lo) & (b2 < hi)
    n = hi - lo
    raw = np.zeros((n, n)); np.maximum.at(raw, (b1[k] - lo, b2[k] - lo), c[k].astype(float))
    raw = np.maximum(raw, raw.T)
    ww = w[lo:hi].astype(float)
    bal = raw * ww[:, None] * ww[None, :]
    bal[np.isnan(ww), :] = np.nan; bal[:, np.isnan(ww)] = np.nan
    return raw.astype(np.float32), bal.astype(np.float32), ww.astype(np.float32), bs[lo:hi].astype(np.int64)


out = {'window': np.array([START, END]), 'res': RES}
for tag, fn in FILES.items():
    with h5py.File(fn, 'r') as f:
        raw, bal, w, s = cut(f[f'resolutions/{RES}'])
    out.update({f'{tag}_raw': raw, f'{tag}_bal': bal, f'{tag}_w': w, f'{tag}_starts': s, f'{tag}_file': fn})
    i = int(np.where(s == 13_225_000)[0][0]); j = i + 29          # 145 kb = 29 bins
    print(f'{tag:6s} {fn}: {len(s)} bins | 13.225 x 13.370 = {raw[i, j]:.0f} read pairs (check: H1 291, HFFc6 456)')
np.savez_compressed('stripe_maps.npz', **out)
print('wrote stripe_maps.npz')
