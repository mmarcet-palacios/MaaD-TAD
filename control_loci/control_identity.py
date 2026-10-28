#!/usr/bin/env python3
"""
control_identity.py - how alike are the homologous copies in PRAMEF and in the two control loci?

For each locus the window is cut into 2 kb bins exactly as in "Homology maps" (Materials and methods):
soft-masked repeats and N removed, bins with < 300 unmasked bases dropped, every pair of bins aligned on
both strands (Smith-Waterman, parasail; match +1, mismatch -3, gap open 5, gap extension 2), score divided by
the shorter unmasked bin x 100. Pairs scoring >= 30 are "homologous" (the paper reports 133 for PCDH and 41
for ZNF_chr19, a check that this script reproduces the map). Each homologous pair is then re-aligned with a
traceback and its percent identity is measured over the aligned columns (matches / (matches + mismatches +
gap columns)); coverage = aligned length of the shorter bin.

Output (in this folder):
  identity_pairs.csv     one row per homologous bin pair (locus, bin_i, bin_j, strand, score_pct, identity, coverage)
  identity_summary.txt   per locus: number of pairs, median identity and its interquartile range, median coverage

Usage (Anaconda (base) prompt, from MaaD-TAD\\control_loci):
  python control_identity.py
Needs numpy and parasail (both already used for the homology maps) and internet access to the UCSC
REST API on the first run (sequences are cached in identity_cache\\ afterwards). Uses all cores but one;
expect a few minutes per locus.
"""
import csv, json, os, sys, time, urllib.request
from multiprocessing import Pool, cpu_count

import numpy as np
import parasail

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "identity_cache")
BIN = 2000
MIN_UNMASKED = 300
THRESH = 30.0
GAP_OPEN, GAP_EXT = 5, 2
MATRIX = parasail.matrix_create("ACGT", 1, -3)

# name, chrom, 0-based start, number of 2 kb bins (windows of the paper; README of this folder)
LOCI = [
    ("PRAMEF", "chr1", 12771320, 325),
    ("PCDH", "chr5", 140700000, 450),
    ("ZNF_chr19", "chr19", 57000000, 500),
]

COMP = str.maketrans("ACGT", "TGCA")


def _ssl_context():
    """Trust the Windows certificate store (truststore) or, failing that, the certifi bundle.
    Plain urllib on some Anaconda installs has no root certificates and fails with CERTIFICATE_VERIFY_FAILED."""
    import ssl
    try:
        import truststore
        return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    except ImportError:
        pass
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


SSL_CTX = _ssl_context()


def rc(s):
    return s.translate(COMP)[::-1]


def fetch(chrom, start, end):
    """Soft-masked hg38 sequence from the UCSC REST API, in 200 kb pieces, cached on disk."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "%s_%d_%d.txt" % (chrom, start, end))
    if os.path.exists(path):
        return open(path).read()
    parts, step = [], 200_000
    for s in range(start, end, step):
        e = min(s + step, end)
        url = ("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=%s;start=%d;end=%d"
               % (chrom, s, e))
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=60, context=SSL_CTX) as r:
                    parts.append(json.load(r)["dna"])
                break
            except Exception as exc:
                if attempt == 4:
                    raise
                print("  retrying UCSC (%s)" % exc, flush=True)
                time.sleep(3 * (attempt + 1))
        time.sleep(0.5)
    seq = "".join(parts)
    if len(seq) != end - start:
        sys.exit("UCSC returned %d bases for %s:%d-%d, expected %d" % (len(seq), chrom, start, end, end - start))
    open(path, "w").write(seq)
    return seq


def unmasked_bins(seq, n):
    """Upper-case A/C/G/T only: soft-masked repeats (lower case) and N are removed."""
    out = []
    for i in range(n):
        b = seq[i * BIN:(i + 1) * BIN]
        u = "".join(c for c in b if c in "ACGT")
        out.append(u if len(u) >= MIN_UNMASKED else None)
    return out


_BINS = None


def _init(bins):
    global _BINS
    _BINS = bins


def _row(i):
    """Score bin i against every later bin on both strands; return the homologous pairs."""
    a = _BINS[i]
    if a is None:
        return []
    prof = parasail.profile_create_16(a, MATRIX)
    hits = []
    for j in range(i + 1, len(_BINS)):
        b = _BINS[j]
        if b is None:
            continue
        brc = rc(b)
        s_f = parasail.sw_striped_profile_16(prof, b, GAP_OPEN, GAP_EXT).score
        s_r = parasail.sw_striped_profile_16(prof, brc, GAP_OPEN, GAP_EXT).score
        best, strand, other = (s_f, "+", b) if s_f >= s_r else (s_r, "-", brc)
        short = min(len(a), len(b))
        pct = 100.0 * best / short
        if pct < THRESH:
            continue
        res = parasail.sw_trace_striped_16(a, other, GAP_OPEN, GAP_EXT, MATRIX)
        cig = res.cigar.decode
        cig = cig.decode() if isinstance(cig, bytes) else cig
        ops, num = [], ""
        for ch in cig:
            if ch.isdigit():
                num += ch
                continue
            ops.append((int(num or 1), ch))
            num = ""
        # parasail pads the CIGAR with the unaligned start of one sequence as a long I or D run:
        # these are not part of the local alignment, so leading and trailing gap runs are dropped
        while ops and ops[0][1] in "IDN":
            ops.pop(0)
        while ops and ops[-1][1] in "IDN":
            ops.pop()
        m = sum(k for k, op in ops if op in "=M")
        x = sum(k for k, op in ops if op == "X")
        g = sum(k for k, op in ops if op in "IDN")
        cols = m + x + g
        ident = 100.0 * m / cols if cols else float("nan")
        cover = 100.0 * min(m + x + g, short) / short
        hits.append((i, j, strand, round(pct, 2), round(ident, 2), round(cover, 1)))
    return hits


def main():
    rows, lines = [], []
    procs = max(1, cpu_count() - 1)
    for name, chrom, start, n in LOCI:
        t0 = time.time()
        print("%s: fetching %s:%d-%d" % (name, chrom, start, start + n * BIN), flush=True)
        seq = fetch(chrom, start, start + n * BIN)
        bins = unmasked_bins(seq, n)
        kept = sum(b is not None for b in bins)
        print("  %d of %d bins kept; aligning on %d cores" % (kept, n, procs), flush=True)
        with Pool(procs, initializer=_init, initargs=(bins,)) as pool:
            hits = [h for r in pool.imap(_row, range(n), chunksize=4) for h in r]
        for h in hits:
            rows.append((name,) + h)
        ident = np.array([h[4] for h in hits])
        cover = np.array([h[5] for h in hits])
        far = np.array([h[4] for h in hits if h[1] - h[0] >= 2])
        if len(ident):
            q1, med, q3 = np.percentile(ident, [25, 50, 75])
            line = ("%-10s pairs %4d (non-adjacent %4d)  identity median %.1f%% (IQR %.1f-%.1f%%), max %.2f%%, "
                    ">=97.9%%: %.0f%% of pairs  coverage median %.0f%%  [%.0f s]"
                    % (name, len(hits), len(far), med, q1, q3, ident.max(),
                       100.0 * np.mean(ident >= 97.9), np.median(cover), time.time() - t0))
        else:
            line = "%-10s no homologous pairs" % name
        print("  " + line, flush=True)
        lines.append(line)
    with open(os.path.join(HERE, "identity_pairs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["locus", "bin_i", "bin_j", "strand", "score_pct", "identity_pct", "coverage_pct"])
        w.writerows(rows)
    with open(os.path.join(HERE, "identity_summary.txt"), "w") as f:
        f.write("Homologous 2 kb bin pairs (score >= 30), identity over aligned columns\n")
        f.write("Paper check: PCDH 133 pairs, ZNF_chr19 41 pairs\n\n")
        f.write("\n".join(lines) + "\n")
    print("\nwrote identity_pairs.csv and identity_summary.txt", flush=True)


if __name__ == "__main__":
    main()
