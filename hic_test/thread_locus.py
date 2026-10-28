"""thread_locus.py -- thread_sequence.py for any locus: put the GRCh38 bases on the all-atom blocks.

Same method as ..\\thread_sequence.py (base swapped by C1' + glycosidic-ring anchors, backbone untouched);
only the locus numbering comes from the command line instead of the PRAMEF constants.

Numbering (the chain spans the whole window, see build_locus.py):
    chain A residue r  <->  genomic position start + r          (base = reference, + strand)
    chain B residue r  <->  genomic position end + 1 - r        (base = complement)
    where start = 0-based window start (bin start), end = window end; FASTA = exactly start+1..end.

Run (called by build_locus.py; by hand from the locus folder):
    python ..\\thread_locus.py --indir blocks --outdir blocks_seq --fasta ..\\PWRN_chr15.fa --start 24130000 --end 24440000
"""
import argparse, os, glob, time
import numpy as np
import gemmi
from multiprocessing import Pool

RES = {'A': 'DA', 'C': 'DC', 'G': 'DG', 'T': 'DT'}
COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
BACKBONE = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
ANCH = {'DA': ("C1'", 'N9', 'C4', 'C8'), 'DG': ("C1'", 'N9', 'C4', 'C8'),
        'DC': ("C1'", 'N1', 'C2', 'C6'), 'DT': ("C1'", 'N1', 'C2', 'C6')}
PAIR = {'DA': 'DT', 'DT': 'DA', 'DG': 'DC', 'DC': 'DG'}


def get_sequence(fasta, n):
    s = ''.join(l.strip() for l in open(fasta) if not l.startswith('>')).upper()
    assert len(s) == n, f'sequence length {len(s)} != window {n}'
    return s


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def xyz(res, name):
    a = res.find_atom(name, '*'); return np.array(a.pos.tolist()) if a else None


def templates(model):
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
    P = np.array([T[new][n] for n in ANCH[new]])
    Q = np.array([xyz(res, n) for n in ANCH[res.name]])
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
    path, outdir, seq, start, end = job
    nres = end - start
    tag = os.path.basename(path)[:-4]; t0 = time.time()
    st = gemmi.read_structure(path); st.setup_entities(); model = st[0]
    T = templates(model)
    rows, changed = [], {'A': [], 'B': []}
    for ch in model:
        if ch.name not in ('A', 'B'): continue
        n = same = diff = nn = 0
        for r in ch:
            g = start + r.seqid.num if ch.name == 'A' else end + 1 - r.seqid.num
            b = seq[g - start - 1]
            if ch.name == 'B': b = COMP.get(b, 'N')
            n += 1
            if b not in RES: nn += 1; continue
            if r.name == RES[b]: same += 1; continue
            swap(r, RES[b], T); diff += 1; changed[ch.name].append(r.seqid.num)
        rows.append((ch.name, n, same, diff, nn))
    A = {r.seqid.num: r for r in model['A']}; B = {r.seqid.num: r for r in model['B']}
    pairs = wc = 0; dists = []
    ch_set = set(changed['A']) | {nres + 1 - x for x in changed['B']}
    for ra, r in A.items():
        p = B.get(nres + 1 - ra)
        if p is None: continue
        pairs += 1; wc += PAIR.get(r.name) == p.name
        if ra in ch_set:
            a1 = xyz(r, 'N1' if r.name in ('DA', 'DG') else 'N3'); a2 = xyz(p, 'N3' if r.name in ('DA', 'DG') else 'N1')
            if a1 is not None and a2 is not None: dists.append(float(np.linalg.norm(a1 - a2)))
    ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
    clash = 0; dmin = 99.0
    for chn in ('A', 'B'):
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
        f.write('chain,residues,already_correct,changed,N\n')
        for row in rows: f.write(','.join(map(str, row)) + '\n')
    tot = sum(r[1] for r in rows); nn = sum(r[4] for r in rows)
    match = sum(r[2] + r[3] for r in rows) / max(1, tot - nn)
    med = np.median(dists) if dists else float('nan')
    return (f'{tag}: {tot} nt, changed {sum(r[3] for r in rows)}, N {nn}, match {match:.4f}, '
            f'WC pairs {wc}/{pairs}, clashes<2.2A {clash} (min {dmin if clash else float("nan"):.2f}), '
            f'median pair N-N {med:.2f} A, {time.time() - t0:.0f}s')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--indir', default='blocks'); ap.add_argument('--outdir', default='blocks_seq')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--start', type=int, required=True, help='0-based window start (bin start)')
    ap.add_argument('--end', type=int, required=True)
    ap.add_argument('--blocks', type=int, nargs='*'); ap.add_argument('--workers', type=int, default=12)
    a = ap.parse_args()
    seq = get_sequence(a.fasta, a.end - a.start)
    print(f'sequence: {len(seq):,} bp, N: {seq.count("N"):,}')
    files = sorted(glob.glob(os.path.join(a.indir, 'block_??.cif')))
    if a.blocks: files = [f for f in files if int(os.path.basename(f)[6:8]) in a.blocks]
    lines = []
    with Pool(a.workers) as p:
        for line in p.imap(do_block, [(f, a.outdir, seq, a.start, a.end) for f in files]):
            print(line, flush=True); lines.append(line)
    with open(os.path.join(a.outdir, 'sequence_summary.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'wrote {len(lines)} blocks to {a.outdir}/')


if __name__ == '__main__':
    main()
