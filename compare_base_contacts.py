"""compare_base_contacts.py -- close contacts of DNA base atoms, before vs after sequence threading.

Run in MaaD-TAD (needs gemmi, numpy):
    python compare_base_contacts.py                 # all blocks, 6 workers
    python compare_base_contacts.py --blocks 30 67  # quick test
    python compare_base_contacts.py --before blocks_seq --after blocks_final   # after tail repair
Compares blocks_tethered/ (template bases) with blocks_seq/ (GRCh38 bases), counted the same way:
every atom pair closer than 2.2 A (0.22 nm) in which at least one atom is a base atom of chain A or B
and the two atoms are in different residues. Pairs are counted once and split by partner:
  DNA-DNA       base atom vs another nucleotide (stacking or pairing neighbours, other duplexes)
  DNA-histone   base atom vs any protein atom
Also reports the number of bases with at least one such contact and the closest distance.
Writes blocks_seq/base_contacts_before_after.csv and prints totals.
"""
import argparse, glob, os, time
from multiprocessing import Pool
import gemmi

BACKBONE = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
CUT = 2.2
BEFORE, AFTER = 'blocks_tethered', 'blocks_seq'


def count(path):
    st = gemmi.read_structure(path); m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, 5).populate()
    seen = set(); dd = dh = 0; bases = set(); dmin = 99.0
    for ci, ch in enumerate(m):
        if ch.name not in ('A', 'B'): continue
        for ri, r in enumerate(ch):
            for ai, a in enumerate(r):
                if a.name in BACKBONE or a.element.name == 'H': continue
                for mk in ns.find_atoms(a.pos, '\0', radius=CUT):
                    if mk.chain_idx == ci and mk.residue_idx == ri: continue
                    cra = mk.to_cra(m)
                    if cra.atom.element.name == 'H': continue
                    d = cra.atom.pos.dist(a.pos)
                    if d >= CUT: continue
                    key = tuple(sorted(((ci, ri, ai), (mk.chain_idx, mk.residue_idx, mk.atom_idx))))
                    if key in seen: continue
                    seen.add(key); bases.add((ci, ri)); dmin = min(dmin, d)
                    if cra.chain.name in ('A', 'B'): dd += 1
                    else: dh += 1
    return dd, dh, len(bases), dmin


def job(arg):
    b, before, after = arg          # passed explicitly: Windows workers do not inherit globals set in main()
    t = time.time(); out = [b]
    for d in (before, after):
        out += list(count(os.path.join(d, f'block_{b:02d}.cif')))
    out.append(round(time.time() - t))
    return out


def main():
    global BEFORE, AFTER
    ap = argparse.ArgumentParser(); ap.add_argument('--blocks', type=int, nargs='*'); ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--before', default='blocks_tethered'); ap.add_argument('--after', default='blocks_seq')
    a = ap.parse_args(); BEFORE, AFTER = a.before, a.after
    bl = a.blocks or sorted(int(os.path.basename(f)[6:8]) for f in glob.glob(os.path.join(AFTER, 'block_??.cif')))
    cols = ['block', 'before_DNA_DNA', 'before_DNA_histone', 'before_bases', 'before_min_A',
            'after_DNA_DNA', 'after_DNA_histone', 'after_bases', 'after_min_A', 'seconds']
    rows = []
    with Pool(a.workers) as p:
        for r in p.imap(job, [(b, BEFORE, AFTER) for b in bl]):
            rows.append(r); print(' '.join(f'{c}={v:.2f}' if isinstance(v, float) else f'{c}={v}' for c, v in zip(cols, r)), flush=True)
    rows.sort()
    out = os.path.join(AFTER, f'base_contacts_{os.path.basename(BEFORE)}_vs_{os.path.basename(AFTER)}.csv')
    with open(out, 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows: f.write(','.join(f'{v:.2f}' if isinstance(v, float) else str(v) for v in r) + '\n')
    S = lambda i: sum(r[i] for r in rows)
    print(f'\nTOTAL over {len(rows)} blocks (contacts < {CUT} A involving a base atom)')
    print(f'  before: DNA-DNA {S(1):,}  DNA-histone {S(2):,}  bases involved {S(3):,}  closest {min(r[4] for r in rows):.2f} A')
    print(f'  after : DNA-DNA {S(5):,}  DNA-histone {S(6):,}  bases involved {S(7):,}  closest {min(r[8] for r in rows):.2f} A')
    print(f'  (before = {BEFORE}, after = {AFTER})')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
