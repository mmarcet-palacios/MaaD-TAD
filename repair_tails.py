"""repair_tails.py -- move histone tails out of the DNA, by torsions only (Builder repair step).

Why: in the built blocks one H3 N-terminal tail per nucleosome (copy E, residues 1-30) runs
through the linker DNA entering the core. The tails are disordered in cells; 1KX5 captured one
conformation, and placing that copy on every core puts it where the crystal had no linker.
How: the same torsion-space relaxation as tail_relax.py. Only backbone phi/psi of the tail rotate,
anchor outward, so bond lengths and bond angles stay exactly as built; the histone fold, the DNA
and all other atoms never move. A tail is moved only if one of its atoms lies within 2.2 A of DNA,
so tails that merely touch DNA (as they do in the crystal) are left alone; the search keeps the
tail away from DNA and from all other protein atoms with equal weight.

Run in MaaD-TAD (needs gemmi, numpy, scipy):
    python repair_tails.py                         # blocks_seq -> blocks_final, 6 workers
    python repair_tails.py --blocks 30 67          # quick test
    python repair_tails.py --indir blocks_final --blocks 7 9 12 --step 5 --rounds 6   # second, finer pass
Writes blocks_final/block_XX.cif and blocks_final/tail_repair.csv, and prints per-block counts of
atom pairs closer than 2.2 A (0.22 nm) between tail atoms and DNA, before and after.
Then rerun the contact count on the repaired blocks:
    python compare_base_contacts.py --after blocks_final
"""
import argparse, glob, os, time
from multiprocessing import Pool
import numpy as np
import gemmi
from scipy.spatial import cKDTree

# tails: chain letter -> list of (side, anchor residue number); anchor and the fold never move
TAILS = {'A': [('N', 40)], 'E': [('N', 40)],            # H3
         'B': [('N', 24)], 'F': [('N', 24)],            # H4
         'C': [('N', 16), ('C', 118)], 'G': [('N', 16), ('C', 118)],   # H2A
         'D': [('N', 32)], 'H': [('N', 32)]}            # H2B
CLASH = 2.2
DNA_TARGET, PROT_TARGET = 3.6, 3.2
CTX_TRIGGER = 4.0          # cross-block mode: a tail within 4 A of a neighbouring block is moved (QC limits: CA-P 4 A, CA-CA 3.8 A)
STEP, ROUNDS = 10.0, 3


def rot(axis, t):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)


def torsions(names, anchor, side, n):
    """(axis_a, axis_b, moving) with atoms as (residue index, name); names[k] = atom names of residue k."""
    out = []
    if side == 'N':
        for j in range(anchor, 0, -1):
            if 'N' in names[j] and 'CA' in names[j]:
                out.append(((j, 'N'), (j, 'CA'), [(k, a) for k in range(j) for a in names[k]]))
            p = names[j - 1]
            if 'CA' in p and 'C' in p:
                mv = [(k, a) for k in range(j - 1) for a in names[k]] + [(j - 1, a) for a in p if a not in ('C', 'O')]
                out.append(((j - 1, 'CA'), (j - 1, 'C'), mv))
    else:
        for j in range(anchor, n - 1):
            if 'CA' in names[j] and 'C' in names[j]:
                out.append(((j, 'CA'), (j, 'C'), [(k, a) for k in range(j + 1, n) for a in names[k]]))
            q = names[j + 1]
            if 'N' in q and 'CA' in q:
                mv = [(k, a) for k in range(j + 2, n) for a in names[k]] + [(j + 1, a) for a in q if a != 'N']
                out.append(((j + 1, 'N'), (j + 1, 'CA'), mv))
    return [t for t in out if t[2]]


def repair_block(job):
    src, dst, step, rounds, ctx_files, dna_t, prot_t = job
    t0 = time.time()
    st = gemmi.read_structure(src); m = st[0]
    # flat arrays of every atom, so trees can be rebuilt from current coordinates
    xyz, is_dna, owner = [], [], []
    for ci, ch in enumerate(m):
        for ri, r in enumerate(ch):
            for a in r:
                xyz.append(a.pos.tolist()); is_dna.append(ch.name in ('A', 'B')); owner.append((ci, ri))
    n_own = len(xyz)
    # context blocks (cross-block mode): their atoms are fixed obstacles, never moved or written
    for cf in ctx_files:
        for ch in gemmi.read_structure(cf)[0]:
            for r in ch:
                for a in r:
                    xyz.append(a.pos.tolist()); is_dna.append(ch.name in ('A', 'B'))
    xyz = np.array(xyz); is_dna = np.array(is_dna)
    ctx_tree = cKDTree(xyz[n_own:]) if len(xyz) > n_own else None
    dna_idx = np.where(is_dna)[0]; dna_tree = cKDTree(xyz[dna_idx])
    prot_idx = np.where(~is_dna)[0]
    start = {}; k = 0
    for ci, ch in enumerate(m):
        for ri, r in enumerate(ch):
            start[(ci, ri)] = k; k += len(r)
    rows = []; fixed = 0
    for ci, ch in enumerate(m):
        if ch.name[:1] != 'P' or ch.name[-1] not in TAILS: continue
        nums = [r.seqid.num for r in ch]; names = [[a.name for a in r] for r in ch]
        for side, anum in TAILS[ch.name[-1]]:
            if anum not in nums: continue
            anc = nums.index(anum)
            tail_res = range(0, anc) if side == 'N' else range(anc + 1, len(ch))
            if not len(tail_res): continue
            gidx = {(ri, nm): start[(ci, ri)] + ai for ri in range(len(ch)) for ai, nm in enumerate(names[ri])}
            tail_atoms = np.array([gidx[(ri, nm)] for ri in tail_res for nm in names[ri]])
            # protein atoms that may not overlap the tail: all others except this tail and +-3 residues of the anchor
            near_anchor = {gidx[(ri, nm)] for ri in range(max(0, anc - 3), min(len(ch), anc + 4)) for nm in names[ri]}
            excl = set(tail_atoms.tolist()) | near_anchor
            pc = xyz[tail_atoms]; ctr = xyz[gidx[(anc, 'CA')]]
            reach = float(np.linalg.norm(pc - ctr, axis=1).max()) + 25.0
            oth = prot_idx[np.linalg.norm(xyz[prot_idx] - ctr, axis=1) < reach]
            oth = oth[~np.isin(oth, np.fromiter(excl, int))]
            loc_dna = dna_idx[np.linalg.norm(xyz[dna_idx] - ctr, axis=1) < reach]
            dtree = cKDTree(xyz[loc_dna]) if len(loc_dna) else None
            ptree = cKDTree(xyz[oth]) if len(oth) else None
            def score_pts(q):
                d = min(float(dtree.query(q)[0].min()), dna_t) if dtree else dna_t
                p = min(float(ptree.query(q)[0].min()), prot_t) if ptree else prot_t
                return min(d, p) + 0.1 * (d + p)      # protein overlaps count as much as DNA ones
            def n_clash(q):
                c = int((dtree.query(q, distance_upper_bound=CLASH)[0] < CLASH).sum()) if dtree else 0
                return c
            d0 = float(dtree.query(pc)[0].min()) if dtree else 99.0
            p0 = float(ptree.query(pc)[0].min()) if ptree else 99.0
            c0 = n_clash(pc)
            near_ctx = ctx_tree is not None and bool((ctx_tree.query(pc, distance_upper_bound=CTX_TRIGGER)[0] < CTX_TRIGGER).any())
            if c0 == 0 and not near_ctx:             # only tails that overlap DNA (< 2.2 A), or come near a context block, move
                continue
            for _ in range(rounds):
                for (ja, na), (jb, nb), mv in torsions(names, anc, side, len(ch)):
                    pivot = xyz[gidx[(jb, nb)]].copy(); axis = pivot - xyz[gidx[(ja, na)]]
                    if np.linalg.norm(axis) < 1e-6: continue
                    ii = np.array([gidx[x] for x in mv]); rel = xyz[ii] - pivot
                    best_t, best_s = 0.0, score_pts(xyz[ii])
                    for t in np.radians(np.arange(-180.0, 180.0, step)):
                        if t == 0: continue
                        s = score_pts(rel @ rot(axis, t).T + pivot)
                        if s > best_s + 1e-9: best_s, best_t = s, t
                    if best_t: xyz[ii] = rel @ rot(axis, best_t).T + pivot
            pc1 = xyz[tail_atoms]
            d1 = float(dtree.query(pc1)[0].min()) if dtree else 99.0
            p1 = float(ptree.query(pc1)[0].min()) if ptree else 99.0
            rows.append((ch.name, side, len(tail_res), c0, n_clash(pc1), round(d0, 2), round(d1, 2), round(p0, 2), round(p1, 2)))
            fixed += 1
    # write back
    k = 0                                            # own atoms only; context atoms are never written
    for ch in m:
        for r in ch:
            for a in r:
                a.pos = gemmi.Position(*xyz[k]); k += 1
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + '.part'; st.make_mmcif_document().write_file(tmp); os.replace(tmp, dst)
    tag = os.path.basename(src)[:-4]
    c0 = sum(r[3] for r in rows); c1 = sum(r[4] for r in rows)
    worst = min([r[6] for r in rows], default=float('nan'))
    line = (f'{tag}: tails relaxed {fixed}, tail-DNA atoms <2.2 A {c0} -> {c1}, '
            f'closest tail-DNA after {worst:.2f} A, {time.time() - t0:.0f}s')
    return tag, rows, line


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', default='blocks_seq'); ap.add_argument('--outdir', default='blocks_final')
    ap.add_argument('--blocks', type=int, nargs='*'); ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--step', type=float, default=STEP, help='torsion grid in degrees (default 10)')
    ap.add_argument('--rounds', type=int, default=ROUNDS, help='sweeps over the tail torsions (default 3)')
    ap.add_argument('--from-qc', help='qc_problems.csv from stage4_qc.py: repair the blocks it names against their partner blocks')
    ap.add_argument('--size', type=int, default=50, help='nucleosomes per block (for --from-qc)')
    a = ap.parse_args()
    if a.from_qc:
        return cross_block(a)
    files = sorted(glob.glob(os.path.join(a.indir, 'block_??.cif')))
    if a.blocks: files = [f for f in files if int(os.path.basename(f)[6:8]) in a.blocks]
    jobs = [(f, os.path.join(a.outdir, os.path.basename(f)), a.step, a.rounds, [], DNA_TARGET, PROT_TARGET) for f in files]
    allrows = []; lines = []
    with Pool(a.workers) as p:
        for tag, rows, line in p.imap(repair_block, jobs):
            print(line, flush=True); lines.append(line); allrows += [(tag,) + r for r in rows]
    os.makedirs(a.outdir, exist_ok=True)
    c0 = sum(r[4] for r in allrows); c1 = sum(r[5] for r in allrows)
    left = sum(1 for r in allrows if r[5] > 0)
    csvp = os.path.join(a.outdir, 'tail_repair.csv')
    rows = [tuple(map(str, r)) for r in sorted(allrows)]
    if a.blocks and os.path.exists(csvp):          # a partial rerun keeps the earlier rows and adds its own as a later pass
        rows = open(csvp).read().splitlines()[1:] + [','.join((r[0] + '_rerun',) + r[1:]) for r in rows]
    else:
        rows = [','.join(r) for r in rows]
    with open(csvp, 'w') as f:
        f.write('block,chain,side,tail_residues,dna_atoms_lt2.2_before,dna_atoms_lt2.2_after,'
                'min_dna_A_before,min_dna_A_after,min_protein_A_before,min_protein_A_after\n')
        f.write('\n'.join(rows) + '\n')
    print(f'\nTOTAL {len(lines)} blocks: {len(allrows)} tails relaxed; tail atoms within 2.2 A of DNA {c0:,} -> {c1:,}; '
          f'tails still touching DNA: {left}')
    print(f'wrote {a.outdir}/ (block_XX.cif, tail_repair.csv)')


def cross_block(a):
    """Tails that the per-block pass swung into a NEIGHBOURING block (found by stage4_qc.py's long-range check).
    Each named block is relaxed again with its partner blocks loaded as fixed obstacles, one block at a time
    (so a partner repaired first is seen in its new form), aiming for 4.5 A from DNA and 4.0 A from protein."""
    import csv
    pairs = set()
    for r in csv.DictReader(open(a.from_qc)):
        bi, bj = int(r['nuc_i']) // a.size, int(r['nuc_j']) // a.size
        if bi != bj: pairs.add((bi, bj)); pairs.add((bj, bi))
    ctx = {}
    for bi, bj in pairs: ctx.setdefault(bi, set()).add(bj)
    blk = lambda b: os.path.join(a.indir, f'block_{b:02d}.cif')
    print('cross-block repair: ' + '; '.join(f'block {b} vs {sorted(c)}' for b, c in sorted(ctx.items())))
    for b in sorted(ctx):
        out = os.path.join(a.outdir, f'block_{b:02d}.cif')
        tag, rows, line = repair_block((blk(b), out, a.step, a.rounds, [blk(c) for c in sorted(ctx[b])], 4.5, 4.0))
        print(line, '| tails moved:', ', '.join(f'{r[0]}{r[1]}' for r in rows) or 'none', flush=True)
    print('done; rerun stage4_qc.py on', a.outdir)


if __name__ == '__main__':
    main()
