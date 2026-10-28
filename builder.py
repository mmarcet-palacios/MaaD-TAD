"""builder.py -- MaaD-TAD module 3, Builder: build the DNA and histones in atoms, close every linker,
thread the GRCh38 sequence and repair the histone tails.

Steps (settings as in S2 Table):
  build     run_stage2.py              full-atom blocks of --size nucleosomes (default 50), built in parallel;
                                       every linker closed to 1e-10 nm -> blocks_tethered/ (about an hour on 24 cores)
  seq       thread_sequence.py         GRCh38 bases onto the blocks (backbone unchanged) -> blocks_seq/
            compare_base_contacts.py   close atom pairs before and after threading (blocks_tethered vs blocks_seq)
  tails     repair_tails.py            histone tails moved by backbone torsions only (10 deg grid, 3 sweeps) -> blocks_final/;
                                       blocks with a tail still touching DNA get a second pass (5 deg, 6 sweeps)
            compare_base_contacts.py   blocks_seq vs blocks_final

Order used for the paper (RUN_PRAMEF.bat): build, then Checker (qc, repair), then seq and tails, then Checker final.

Run from the MaaD-TAD folder:
    python builder.py --steps build
    python builder.py --steps seq tails
    python builder.py --steps build --size 25      # smaller blocks (block size is a practical choice, not a structural one)
"""
import csv, os
from maadtad_steps import ROOT, run, main, workers

FRAMES = 'frames_tethered.npz'


def build(a):
    run('run_stage2.py', '--frames', FRAMES, '--outdir', 'blocks_tethered', '--size', a.size,
        '--workers', workers(a, 24), dry=a.dry_run)


def seq(a):
    run('thread_sequence.py', '--indir', 'blocks_tethered', '--outdir', 'blocks_seq', '--fasta', 'locus_GRCh38.fa',
        '--workers', workers(a, 6), dry=a.dry_run)
    run('compare_base_contacts.py', '--before', 'blocks_tethered', '--after', 'blocks_seq',
        '--workers', workers(a, 6), dry=a.dry_run)


def still_touching(path):
    if not os.path.exists(path):
        return []
    rows = csv.DictReader(open(path))
    return sorted({int(r['block'][6:8]) for r in rows if int(r['dna_atoms_lt2.2_after']) > 0})


def tails(a):
    run('repair_tails.py', '--indir', 'blocks_seq', '--outdir', 'blocks_final', '--size', a.size,
        '--workers', workers(a, 6), dry=a.dry_run)
    left = [] if a.dry_run else still_touching(os.path.join(ROOT, 'blocks_final', 'tail_repair.csv'))
    if left:
        print(f'  blocks with a tail still within 0.22 nm of DNA: {left} -> finer pass')
        run('repair_tails.py', '--indir', 'blocks_final', '--outdir', 'blocks_final', '--blocks', *left,
            '--step', 5, '--rounds', 6, '--workers', workers(a, 6), dry=a.dry_run)
    run('compare_base_contacts.py', '--before', 'blocks_seq', '--after', 'blocks_final',
        '--workers', workers(a, 6), dry=a.dry_run)


def options(ap):
    ap.add_argument('--size', type=int, default=50, help='nucleosomes per block (default 50, as in the paper)')


if __name__ == '__main__':
    main('builder', __doc__, {'build': build, 'seq': seq, 'tails': tails}, options)
