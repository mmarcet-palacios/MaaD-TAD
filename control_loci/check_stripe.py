"""check_stripe.py -- is the 145 kb parallel line in Fig 1B present in other Micro-C / Hi-C maps?

Reads one or more .mcool files (h5py only) at 5 kb and reports, for the D3/D4 line
(rows 13.225-13.245 Mb against columns +145 kb), the raw read pairs on the line,
the mean of the bins one step off the line, and the local (diagonal) contacts of each bin.

  python check_stripe.py H1_4DNFI9GMP2J8.mcool HFFc6_microC.mcool  [more .mcool files]

Same coordinates as Fig 1B (GRCh38).
"""
import sys, h5py, numpy as np
CHROM, RES, OFF = "chr1", 5000, 145000
ROWS = list(range(13215000, 13255000, 5000))           # 13.215-13.250 Mb (line runs 13.225-13.245)
def block(path, lo, hi):
    f = h5py.File(path, "r"); g = f[f"resolutions/{RES}"] if "resolutions" in f else f
    names = [c.decode() if isinstance(c, bytes) else c for c in g["chroms/name"][:]]
    cname = next(n for n in names if n in (CHROM, CHROM.replace("chr", ""))); cid = names.index(cname)
    bc, bs = g["bins/chrom"][:], g["bins/start"][:]
    w = g["bins/weight"][:] if "weight" in g["bins"] else np.ones(len(bs))
    sel = np.where((bc == cid) & (bs >= lo) & (bs < hi))[0]; a, b = sel[0], sel[-1] + 1
    off = g["indexes/bin1_offset"]; p0, p1 = off[a], off[b]
    b1, b2, c = g["pixels/bin1_id"][p0:p1], g["pixels/bin2_id"][p0:p1], g["pixels/count"][p0:p1]
    k = (b2 >= a) & (b2 < b); M = np.zeros((b - a, b - a))
    np.maximum.at(M, (b1[k] - a, b2[k] - a), c[k]); M = np.maximum(M, M.T)
    return M, bs[a:b], w[a:b]
for path in sys.argv[1:]:
    M, st, w = block(path, 13200000, 13420000)
    ix = {int(s): i for i, s in enumerate(st)}; tot = 0
    print(f"\n{path}\n  row(Mb)  col(Mb)   on-line  off-line(mean of 4)  local_row  local_col  weight_row  weight_col")
    for r in ROWS:
        i, j = ix[r], ix[r + OFF]
        offl = np.mean([M[i, j - 1], M[i, j + 1], M[i - 1, j], M[i + 1, j]])
        loc = lambda k: M[k, k] + M[k, k + 1] + M[k - 1, k]
        print(f"  {r/1e6:7.3f}  {(r+OFF)/1e6:7.3f}  {M[i,j]:8.0f}  {offl:19.1f}  {loc(i):9.0f}  {loc(j):9.0f}  {w[i]:10.4g}  {w[j]:10.4g}")
        tot += M[i, j]
    print(f"  total on line (13.215-13.250): {tot:.0f} read pairs")
