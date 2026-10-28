#!/usr/bin/env python3
"""
misassign_pramef.py - can reads assigned to the wrong copy produce the 145 kb PRAMEF line?

The line joins row bins chr1:13.225-13.245 Mb to column bins 145 kb downstream (13.370-13.390 Mb).
A read pair from a short-range contact inside one copy lands on the line only if one mate is placed in the
other copy. For every 5 kb bin of the line (both copies) this script
  1. aligns the bin to its partner copy (minimap2, asm20) and reports the identity and substitution count;
  2. simulates reads of 50, 100 and 150 nt every 10 nt along the bin, with 0.5% random substitutions, maps them
     to the whole PRAMEF locus (minimap2 short-read preset) and counts reads whose best hit falls in the partner
     bin at MAPQ >= 1 and >= 30 (per-read misassignment rate m);
  3. repeats step 2 on a "converted" bin in which a fraction p of the substitutions that separate the two copies
     carry the partner's base (a model of gene conversion or a reference error in this individual);
  4. turns m into an expected number of read pairs on each line cell:
        expected = 2 x (m_row x local_row + m_col x local_col)
     where local_* is the raw count of read pairs inside the bin (the diagonal cell; the short-range contacts that
     would be moved onto the line) in the H1 Micro-C map, and the factor 2 counts either mate being misplaced.
     For a conversion carried on one of the two homologs, local counts are halved.
Inputs : locus_GRCh38.fa (chr1:12,771,346-13,423,345), control_loci/fig1b_microc_H1.npz (H1 raw counts at 5 kb)
Output : misassign_report.txt, misassign_table.csv
Usage (from MaaD-TAD/control_loci; needs numpy, mappy):
  python misassign_pramef.py ../locus_GRCh38.fa fig1b_microc_H1.npz
"""
import sys, random, csv
import numpy as np
import mappy as mp

FA = sys.argv[1] if len(sys.argv) > 1 else 'locus_GRCh38.fa'
NPZ = sys.argv[2] if len(sys.argv) > 2 else 'fig1b_microc_H1.npz'
L0 = 12_771_346                       # 1-based start of the locus FASTA
ROWS = [13_225_000, 13_230_000, 13_235_000, 13_240_000, 13_245_000]
OFF = 145_000
LENS = (50, 100, 150)
STEP = 10
ERR = 0.005
PS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
OBS_H1 = {13_225_000: 291, 13_230_000: 197, 13_235_000: 129, 13_240_000: 122, 13_245_000: 55}
rng = random.Random(1)

seq = ''.join(l.strip() for l in open(FA) if not l.startswith('>')).upper()
sub = lambda a, b: seq[a - L0:b - L0]
ref = mp.Aligner(seq=seq, preset='sr')


def partner_alignment(a, b, pa, pb):
    """align bin [a,b) to partner window [pa,pb); return identity, list of (q_pos, partner_base) substitutions"""
    al = mp.Aligner(seq=sub(pa, pb), preset='asm20')
    best = None
    for h in al.map(sub(a, b), cs=True):
        if h.strand == 1 and (best is None or h.mlen > best.mlen):
            best = h
    if best is None:
        return None, []
    subs, q = [], best.q_st
    import re
    for op in re.findall(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+)', best.cs):
        if op[0] == ':':
            q += int(op[1:])
        elif op[0] == '*':
            subs.append((q, op[1].upper())); q += 1    # *ref query: op[1] = partner (reference) base
        elif op[0] == '+':
            q += len(op) - 1
    return best.mlen / best.blen, subs


def mutate(s):
    s = list(s)
    for i in range(len(s)):
        if rng.random() < ERR:
            s[i] = rng.choice([c for c in 'ACGT' if c != s[i]])
    return ''.join(s)


def misrate(src, partner_lo, partner_hi, L):
    n = n1 = n30 = 0
    for i in range(0, len(src) - L + 1, STEP):
        r = src[i:i + L]
        if 'N' in r:
            continue
        n += 1
        hits = list(ref.map(mutate(r)))
        if not hits:
            continue
        h = hits[0]
        pos = h.r_st + L0
        if partner_lo <= pos < partner_hi:
            if h.mapq >= 1:
                n1 += 1
            if h.mapq >= 30:
                n30 += 1
    return n, n1 / max(n, 1), n30 / max(n, 1)


d = np.load(NPZ)
starts, raw = d['starts_5000'], d['raw_5000']
idx = {int(s): k for k, s in enumerate(starts)}

rows_out, lines = [], []
lines.append('MISASSIGNMENT TEST, 145 kb PRAMEF line (H1 Micro-C, GRCh38)')
lines.append(f'reads every {STEP} nt, lengths {LENS}, {ERR*100:.1f}% substitutions; mapped to the whole locus (minimap2 sr)')
for r in ROWS:
    c = r + OFF
    for (a, pa, label) in ((r, c, 'row copy -> column copy'), (c, r, 'column copy -> row copy')):
        ident, subs = partner_alignment(a, a + 5000, pa - 2000, pa + 7000)
        for p in PS:
            src = list(sub(a, a + 5000))
            k = int(round(p * len(subs)))
            for q, base in rng.sample(subs, k):
                if q < len(src):
                    src[q] = base
            src = ''.join(src)
            for L in LENS:
                n, m1, m30 = misrate(src, pa - 1000, pa + 6000, L)
                rows_out.append(dict(row_bin=r, source_bin=a, direction=label, identity=round(ident or 0, 4),
                                     n_subs=len(subs), converted_fraction=p, read_len=L, reads=n,
                                     mis_mapq1=round(m1, 5), mis_mapq30=round(m30, 5)))
    print('done', r, flush=True)

with open('misassign_table.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)

# expected read pairs on the line
lines.append('')
lines.append('per bin: identity to partner, substitutions, local read pairs (H1 raw, diagonal cell), observed on line')
loc = {}
for r in ROWS:
    c = r + OFF
    lr = float(raw[idx[r], idx[r]]); lc = float(raw[idx[c], idx[c]])
    loc[r] = (lr, lc)
    idr = next(x for x in rows_out if x['source_bin'] == r)
    idc = next(x for x in rows_out if x['source_bin'] == c)
    lines.append(f'  {r/1e6:.3f} Mb  identity {idr["identity"]*100:.1f}% ({idr["n_subs"]} subs) | local {lr:.0f} + {lc:.0f} | observed {OBS_H1[r]}')
lines.append('')
lines.append('expected read pairs on the line (sum over the 5 cells) by converted fraction p, read length and MAPQ filter')
lines.append('   p    len   MAPQ>=1 (both homologs)  MAPQ>=30 (both)  MAPQ>=30 (one homolog)   observed')
obs = sum(OBS_H1.values())
for p in PS:
    for L in LENS:
        e1 = e30 = 0.0
        for r in ROWS:
            c = r + OFF
            mr = next(x for x in rows_out if x['source_bin'] == r and x['converted_fraction'] == p and x['read_len'] == L)
            mc = next(x for x in rows_out if x['source_bin'] == c and x['converted_fraction'] == p and x['read_len'] == L)
            lr, lc = loc[r]
            e1 += 2 * (mr['mis_mapq1'] * lr + mc['mis_mapq1'] * lc)
            e30 += 2 * (mr['mis_mapq30'] * lr + mc['mis_mapq30'] * lc)
        lines.append(f'  {p:4.2f}  {L:4d}   {e1:12.0f}            {e30:10.0f}        {e30/2:10.0f}            {obs}')
open('misassign_report.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
