"""thread_sequence.py -- Builder step: put the real GRCh38 sequence on the all-atom model.

The blocks in blocks_tethered/ carry the right backbone path and residue numbering, but their
bases come from the build templates (1KX5 alpha-satellite on the cores, template bases on the
linkers). This script replaces every base with the GRCh38 base at that position, keeping the
sugar-phosphate backbone exactly where it is.

Numbering (same as stage2_block.py):
    chain A residue r  <->  genomic position r + 12,771,345   (base = reference, + strand)
    chain B residue r  <->  genomic position 13,423,346 - r   (base = complement)

How a base is swapped: the new base (template from the model itself) is superposed on the old
one by C1' and the three glycosidic-ring anchors (purine N9, C4, C8 / pyrimidine N1, C2, C6),
the classic 3DNA-style mutation. The backbone atoms are not touched, so bonds, linker closure
and the block seams are unchanged. Positions inside the 50 kb assembly gap (N in GRCh38)
keep their template base and are counted as 'N'.

Run in MaaD-TAD (needs gemmi, numpy; requests only if the sequence is fetched):
    python thread_sequence.py                     # all blocks, 6 workers
    python thread_sequence.py --blocks 0 1 --workers 2      # quick test
    python thread_sequence.py --fasta locus_GRCh38.fa       # use a local sequence instead of UCSC

Inputs : blocks_tethered/block_XX.cif
Outputs: blocks_seq/block_XX.cif            (new folder; originals untouched)
         blocks_seq/block_XX_seq.csv        per-chain counts and checks
         blocks_seq/sequence_summary.txt    one line per block + totals
         locus_GRCh38.fa                    the locus sequence (cached after the first fetch)
Checks per block:
  - fraction of residues whose base now matches GRCh38 (should be 1.000 outside the gap)
  - Watson-Crick complementarity of every A/B pair after threading
  - new base atoms closer than 2.2 A to any atom outside their own residue (clash count)
  - median N1/N3 pair distance for changed pairs (about 2.9 A for a good pair)
"""
import argparse, os, sys, glob, time
import numpy as np
import gemmi
from multiprocessing import Pool

LOCUS0, LOCUS1 = 12_771_346, 13_423_345          # 1-based inclusive, GRCh38 chr1
OFF_A, OFF_B = 12_771_345, 13_423_346
GAP = (12_954_385, 13_004_384)
RES = {'A': 'DA', 'C': 'DC', 'G': 'DG', 'T': 'DT'}
COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
BACKBONE = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
ANCH = {'DA': ("C1'", 'N9', 'C4', 'C8'), 'DG': ("C1'", 'N9', 'C4', 'C8'),
        'DC': ("C1'", 'N1', 'C2', 'C6'), 'DT': ("C1'", 'N1', 'C2', 'C6')}


def get_sequence(fasta):
    if fasta and os.path.exists(fasta):
        s = ''.join(l.strip() for l in open(fasta) if not l.startswith('>'))
    else:
        import requests
        url = 'https://api.genome.ucsc.edu/getData/sequence'
        r = requests.get(url, params={'genome': 'hg38', 'chrom': 'chr1', 'start': LOCUS0 - 1, 'end': LOCUS1}, timeout=300)
        r.raise_for_status(); s = r.json()['dna']
        with open('locus_GRCh38.fa', 'w') as f:
            f.write(f'>chr1:{LOCUS0}-{LOCUS1} GRCh38 (UCSC hg38)\n')
            for i in range(0, len(s), 80): f.write(s[i:i + 80] + '\n')
    s = s.upper()
    assert len(s) == LOCUS1 - LOCUS0 + 1, f'sequence length {len(s)} != {LOCUS1 - LOCUS0 + 1}'
    return s


def kabsch(P, Q):
    """R, t with R @ P_i + t ~ Q_i."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc); U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def xyz(res, name):
    a = res.find_atom(name, '*'); return np.array(a.pos.tolist()) if a else None


def templates(model):
    """One base of each type, taken from the model itself (anchors + base atoms)."""
    out = {}
    for ch in model:
        if ch.name not in ('A', 'B'): continue
        for r in ch:
            if r.name in ANCH and r.name not in out and all(r.find_atom(n, '*') for n in ANCH[r.name]):
                out[r.name] = {a.name: np.array(a.pos.tolist()) for a in r if a.name not in BACKBONE or a.name == "C1'"}
        if len(out) == 4: break
    assert len(out) == 4, 'could not find all four base types'
    return out


def swap(res, new, T):
    old = res.name
    P = np.array([T[new][n] for n in ANCH[new]])
    Q = np.array([xyz(res, n) for n in ANCH[old]])
    R, t = kabsch(P, Q)
    for nm in [a.name for a in res if a.name not in BACKBONE]:
        res.remove_atom(nm, '*', gemmi.Element(nm[0]))
    for nm, v in T[new].items():
        if nm == "C1'": continue
        at = gemmi.Atom(); at.name = nm; at.element = gemmi.Element(nm[0])
        at.pos = gemmi.Position(*(R @ v + t)); at.occ = 1.0; at.b_iso = 30.0
        res.add_atom(at)
    res.name = new


def do_block(job):
    path, outdir, seq = job
    tag = os.path.basename(path)[:-4]
    t0 = time.time()
    st = gemmi.read_structure(path); st.setup_entities(); model = st[0]
    T = templates(model)
    rows, changed = [], {'A': [], 'B': []}
    for ch in model:
        if ch.name not in ('A', 'B'): continue
        n = same = diff = gapn = 0
        for r in ch:
            g = r.seqid.num + OFF_A if ch.name == 'A' else OFF_B - r.seqid.num
            b = seq[g - LOCUS0]
            if ch.name == 'B': b = COMP.get(b, 'N')
            n += 1
            if b not in RES: gapn += 1; continue
            if r.name == RES[b]: same += 1; continue
            swap(r, RES[b], T); diff += 1; changed[ch.name].append(r.seqid.num)
        rows.append((ch.name, n, same, diff, gapn))
    # checks
    A = {r.seqid.num: r for r in model['A']}; B = {r.seqid.num: r for r in model['B']}
    pairs = wc = 0; dists = []
    ch_set = set(changed['A']) | {652001 - x for x in changed['B']}
    pair_of = {'DA': 'DT', 'DT': 'DA', 'DG': 'DC', 'DC': 'DG'}
    for ra, r in A.items():
        p = B.get(652001 - ra)
        if p is None: continue
        pairs += 1; wc += pair_of.get(r.name) == p.name
        if ra in ch_set:
            a1 = xyz(r, 'N1' if r.name in ('DA', 'DG') else 'N3'); a2 = xyz(p, 'N3' if r.name in ('DA', 'DG') else 'N1')
            if a1 is not None and a2 is not None: dists.append(float(np.linalg.norm(a1 - a2)))
    ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    clash = 0; dmin = 99.0
    for chn in ('A', 'B'):
        ch = model[chn]
        for num in changed[chn]:
            r = A[num] if chn == 'A' else B[num]
            for a in r:
                if a.name in BACKBONE: continue
                for m in ns.find_atoms(a.pos, '\0', radius=2.2):
                    cra = m.to_cra(model)
                    if cra.residue.seqid.num == r.seqid.num and cra.chain.name == chn: continue
                    d = cra.atom.pos.dist(a.pos)
                    if d < 2.2: clash += 1; dmin = min(dmin, d)
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, tag + '.cif'); tmp = out + '.part'
    st.make_mmcif_document().write_file(tmp); os.replace(tmp, out)
    with open(os.path.join(outdir, tag + '_seq.csv'), 'w') as f:
        f.write('chain,residues,already_correct,changed,gap_N\n')
        for row in rows: f.write(','.join(map(str, row)) + '\n')
    tot = sum(r[1] for r in rows); gapn = sum(r[4] for r in rows)
    match = (sum(r[2] + r[3] for r in rows)) / max(1, tot - gapn)
    med = np.median(dists) if dists else float('nan')
    return (f'{tag}: {tot} nt, changed {sum(r[3] for r in rows)}, gap N {gapn}, match {match:.4f}, '
            f'WC pairs {wc}/{pairs}, clashes<2.2A {clash} (min {dmin if clash else float("nan"):.2f}), '
            f'median pair N-N {med:.2f} A, {time.time() - t0:.0f}s')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', default='blocks_tethered'); ap.add_argument('--outdir', default='blocks_seq')
    ap.add_argument('--blocks', type=int, nargs='*'); ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--fasta', default='locus_GRCh38.fa')
    a = ap.parse_args()
    seq = get_sequence(a.fasta)
    print(f'sequence: {len(seq):,} bp, N in gap: {seq.count("N"):,}')
    files = sorted(glob.glob(os.path.join(a.indir, 'block_??.cif')))
    if a.blocks: files = [f for f in files if int(os.path.basename(f)[6:8]) in a.blocks]
    jobs = [(f, a.outdir, seq) for f in files]
    lines = []
    with Pool(a.workers) as p:
        for line in p.imap(do_block, jobs):
            print(line, flush=True); lines.append(line)
    with open(os.path.join(a.outdir, 'sequence_summary.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'wrote {len(lines)} blocks to {a.outdir}/')


if __name__ == '__main__':
    main()
