"""extract_hic_test.py -- inputs for the Hi-C-only test of MaaD-TAD on two gap-free loci (GRCh38).

Loci (2 kb bins, aligned to the Micro-C bin grid so SEQSIM and Micro-C bins coincide exactly):
  PWRN_chr15  chr15:24,130,000-24,440,000  (155 bins; PWRN1/PWRN2 lncRNAs, 7 parallel lines)
  IGHC_chr14  chr14:105,540,000-105,750,000 (105 bins; IGH constant region, 4 parallel lines)

For each locus it writes, in this folder (MaaD-TAD\\hic_test):
  <name>_microc.npz          raw counts, balanced matrix, weights and bin starts, H1 and HFFc6, 2 kb and 5 kb
  <name>_H1_microc.csv       H1 balanced 2 kb matrix in the CSynth layout (Name,Section_1..; same as the SEQSIM csv),
                             scaled so the 99.9th percentile of off-diagonal values = 100, diagonal 100,
                             discarded bins (NaN weight) = 0
  <name>_HFFc6_microc.csv    same for HFFc6
  <name>.fa                  the window sequence (header >chrN:start-end, 1-based)
  <name>_seqsim_{bins.tsv,ff.csv,fr.csv}   raw output of seqsim_windows_fast.exe
  <name>_SEQSIM.csv          SEQSIM map in the CSynth layout: max(direct, inverted), N bins 0.4, diagonal 100
  <name>_preview.png         H1 Micro-C above the diagonal, SEQSIM below (check against the candidate-line figure)
  summary.txt                bins, discarded bins, coverage, N bases

Run (chromatin env, needs numpy, h5py, matplotlib), from MaaD-TAD\\hic_test:
    python extract_hic_test.py
Options: --locus PWRN_chr15 (one locus only), --skip-seqsim (Micro-C only).
Reading the genome FASTA takes a few minutes; SEQSIM takes about a minute per locus.
Reads the .mcool with h5py only; duplicate pixels keep the max, never summed (cooler 0.8.3 quirk).
"""
import argparse, gzip, os, re, subprocess, sys
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
CL = os.path.join(HERE, '..', 'control_loci')
MCOOL = {'H1': os.path.join(CL, '4DNFI9GMP2J8.mcool'), 'HFFc6': os.path.join(CL, '4DNFI643OYP9.mcool')}
FASTA = os.path.join(HERE, '..', 'data', 'GCF_000001405.40_GRCh38.p14_genomic.fna.gz')
SEQSIM_EXE = os.path.join(CL, 'seqsim_windows_fast.exe')
LOCI = {'PWRN_chr15': ('chr15', 24_130_000, 24_440_000),
        'IGHC_chr14': ('chr14', 105_540_000, 105_750_000)}
BIN = 2000


def cut(grp, chrom, start, end):
    names = [c.decode() if isinstance(c, bytes) else c for c in grp['chroms/name'][:]]
    cname = next(n for n in names if n in (chrom, chrom[3:]))
    cid = names.index(cname)
    bc, bs = grp['bins/chrom'][:], grp['bins/start'][:]
    w = grp['bins/weight'][:] if 'weight' in grp['bins'] else np.ones(len(bs))
    sel = np.where((bc == cid) & (bs >= start) & (bs < end))[0]
    lo, hi = sel[0], sel[-1] + 1
    n = hi - lo
    off = grp['indexes/bin1_offset']
    p0, p1 = int(off[lo]), int(off[hi])
    b1 = grp['pixels/bin1_id'][p0:p1]; b2 = grp['pixels/bin2_id'][p0:p1]; c = grp['pixels/count'][p0:p1]
    k = (b2 >= lo) & (b2 < hi)
    b1, b2, c = b1[k] - lo, b2[k] - lo, c[k].astype(np.float64)
    raw = np.zeros((n, n)); np.maximum.at(raw, (b1, b2), c); raw = np.maximum(raw, raw.T)
    ww = w[lo:hi].astype(float)
    bal = raw * ww[:, None] * ww[None, :]
    bal[np.isnan(ww), :] = np.nan; bal[:, np.isnan(ww)] = np.nan
    return raw.astype(np.float32), bal.astype(np.float32), ww.astype(np.float32), bs[lo:hi].astype(np.int64)


def write_csynth(path, M):
    n = len(M)
    with open(path, 'w') as f:
        f.write('Name,' + ','.join(f'Section_{i + 1}' for i in range(n)) + '\n')
        for i in range(n):
            f.write(f'Section_{i + 1},' + ','.join(f'{v:.5f}' for v in M[i]) + '\n')


def microc_to_csynth(bal):
    M = np.nan_to_num(bal.astype(float), nan=0.0)
    iu = np.triu_indices(len(M), 1)
    vals = M[iu][M[iu] > 0]
    s = 100.0 / np.percentile(vals, 99.9) if len(vals) else 1.0
    M = np.minimum(M * s, 100.0)
    np.fill_diagonal(M, 100.0)
    return M


def read_fasta(chroms):
    want = {c[3:] for c in chroms}
    seqs, cur = {}, None
    with gzip.open(FASTA, 'rt') as f:
        for ln in f:
            if ln[0] == '>':
                if cur:
                    seqs[cur] = ''.join(seqs[cur]).upper()
                    if len(seqs) == len(want) and all(isinstance(v, str) for v in seqs.values()):
                        return seqs
                m = re.match(r'>NC_0000\d\d\.\d+ Homo sapiens chromosome (\w+), GRCh38\.p14 Primary Assembly', ln)
                cur = f'chr{m[1]}' if m and m[1] in want else None
                if cur:
                    seqs[cur] = []; print('  reading', cur, flush=True)
            elif cur:
                seqs[cur].append(ln.strip())
    if cur:
        seqs[cur] = ''.join(seqs[cur]).upper()
    return seqs


def preview(path, name, h1, sim, chrom, start, end):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(h1)
    A = np.full((n, n), np.nan)
    iu = np.triu_indices(n, 1); il = np.tril_indices(n, -1)
    lh = np.log10(np.where(h1 > 0, h1, np.nan))
    A[iu] = lh[iu]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(np.where(np.isnan(A), np.nan, A), cmap='Reds', interpolation='nearest',
              vmin=np.nanpercentile(lh[iu], 5), vmax=np.nanpercentile(lh[iu], 99.5))
    B = np.full((n, n), np.nan); ls = np.log10(np.maximum(sim, 0.4)); B[il] = ls[il]
    ax.imshow(B, cmap='Oranges', interpolation='nearest', vmin=np.log10(0.4), vmax=np.log10(20))
    disc = np.where(np.all(np.isnan(h1), 1))[0]
    for d in disc:
        ax.axhline(d, color='#cfd6e0', lw=1); ax.axvline(d, color='#cfd6e0', lw=1)
    t = np.linspace(0, n - 1, 5)
    ax.set_xticks(t); ax.set_xticklabels([f'{(start + i * BIN) / 1e6:.2f}' for i in t])
    ax.set_yticks(t); ax.set_yticklabels([f'{(start + i * BIN) / 1e6:.2f}' for i in t])
    ax.set_title(f'{name}  {chrom}:{start:,}-{end:,}\nH1 Micro-C (above, balanced, log)  ·  SEQSIM (below, log)', fontsize=9)
    ax.set_xlabel(f'{chrom} (Mb)')
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--locus', choices=list(LOCI))
    ap.add_argument('--skip-seqsim', action='store_true')
    a = ap.parse_args()
    loci = {a.locus: LOCI[a.locus]} if a.locus else LOCI
    os.chdir(HERE)
    log = open('summary.txt', 'a')
    def say(*s):
        msg = ' '.join(str(x) for x in s); print(msg, flush=True); log.write(msg + '\n')

    files = {k: h5py.File(v, 'r') for k, v in MCOOL.items()}
    h1_2kb = {}
    for name, (chrom, start, end) in loci.items():
        say(f'== {name}  {chrom}:{start:,}-{end:,}  ({(end - start) // BIN} bins of 2 kb)')
        out = {'chrom': chrom, 'window': np.array([start, end])}
        for cell, f in files.items():
            for res in (2000, 5000):
                raw, bal, w, starts = cut(f[f'resolutions/{res}'], chrom, start, end)
                out.update({f'{cell}_raw_{res}': raw, f'{cell}_bal_{res}': bal,
                            f'{cell}_w_{res}': w, f'{cell}_starts_{res}': starts})
                cov = raw.sum(1)
                say(f'  {cell:5s} {res // 1000} kb: {len(starts)} bins, discarded (NaN weight) {int(np.isnan(w).sum())}, '
                    f'median read pairs per bin {np.median(cov):.0f}')
                if res == 2000:
                    if starts[0] != start or len(starts) != (end - start) // BIN:
                        say('  WARNING: bin grid does not match the window')
                    write_csynth(f'{name}_{cell}_microc.csv', microc_to_csynth(bal))
                    if cell == 'H1':
                        h1_2kb[name] = bal
        np.savez_compressed(f'{name}_microc.npz', **out)
        say(f'  wrote {name}_microc.npz, {name}_H1_microc.csv, {name}_HFFc6_microc.csv')

    if a.skip_seqsim:
        return
    if not os.path.exists(SEQSIM_EXE):
        sys.exit(f'missing {SEQSIM_EXE}: build it with control_loci\\build_seqsim_fast.bat')
    print('reading genome FASTA (a few minutes)...', flush=True)
    seqs = read_fasta({c for c, _, _ in loci.values()})
    for name, (chrom, start, end) in loci.items():
        sub = seqs[chrom][start:end]                      # 0-based start = bin grid of the mcool
        with open(f'{name}.fa', 'w') as o:
            o.write(f'>{chrom}:{start + 1}-{end}\n')
            for i in range(0, len(sub), 80):
                o.write(sub[i:i + 80] + '\n')
        say(f'== {name}: {len(sub):,} bp, {sub.count("N")} N')
        subprocess.run([SEQSIM_EXE, f'{name}.fa', f'{name}_seqsim', str(BIN)], check=True)
        bins = np.genfromtxt(f'{name}_seqsim_bins.tsv', delimiter='\t', skip_header=1, dtype=np.int64).reshape(-1, 3)
        ff = np.loadtxt(f'{name}_seqsim_ff.csv', delimiter=','); fr = np.loadtxt(f'{name}_seqsim_fr.csv', delimiter=',')
        n = (end - start) // BIN
        M = np.full((n, n), 0.4)
        bidx = bins[:, 1]
        M[np.ix_(bidx, bidx)] = np.maximum(np.maximum(ff, fr), 0.4)
        np.fill_diagonal(M, 100.0)
        write_csynth(f'{name}_SEQSIM.csv', M)
        off = M[np.triu_indices(n, 2)]
        say(f'  SEQSIM: {len(bidx)}/{n} bins usable; median {np.median(off):.2f}; pairs > 1.1: {(off > 1.1).sum()}; '
            f'wrote {name}_SEQSIM.csv')
        if name in h1_2kb:
            preview(f'{name}_preview.png', name, h1_2kb[name], M, chrom, start, end)
            say(f'  wrote {name}_preview.png')
    log.close()


if __name__ == '__main__':
    main()
