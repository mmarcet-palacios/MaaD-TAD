"""checker.py -- MaaD-TAD module 4, Checker: run the structural checks on the whole model and rebuild
the blocks that fail them.

Steps (settings as in S2 Table):
  qc       stage4_qc.py on blocks_tethered/: completeness, covalent continuity across all seams,
           phosphodiester geometry (O3'-P-O5' outside 95-115 deg) and long-range clashes between parts more
           than 5 nucleosomes apart (P-P < 0.5 nm, CA-P < 0.4 nm, CA-CA < 0.38 nm) -> qc_summary.txt, qc_problems.csv
  repair   (only if qc_problems.csv lists problems) up to 3 passes of
             relax_frames.py --problems qc_problems.csv   targeted relaxation of the flagged nucleosomes
             tethers_layer3.py                           ring waypoints for the adjusted frames
             check_standin.py                            DNA stand-in against the new waypoints (exit 1 = repeat)
           then pairing_scorecard.py, run_stage2.py --adopt (rebuilds ONLY the blocks whose frames changed)
           and stage4_qc.py again
  final    stage4_qc.py on blocks_final/ (after Builder's seq and tails steps); tails that the per-block
           repair swung into a neighbouring block are relaxed again with the neighbours held fixed
           (repair_tails.py --from-qc), then the check is repeated

Run from the MaaD-TAD folder:
    python checker.py --steps qc repair
    python checker.py --steps final
"""
import os, shutil
from maadtad_steps import ROOT, run, main, problems, workers

FRAMES = 'frames_tethered.npz'


def qc(a, blocks='blocks_tethered'):
    run('stage4_qc.py', '--frames', FRAMES, '--blocks', blocks, '--size', a.size, '--workers', workers(a, 16),
        dry=a.dry_run)


def repair(a):
    if not a.dry_run and problems() == 0:
        print('  qc_problems.csv is empty: nothing to repair')
        return
    # remember what the existing blocks were built from, so that only changed blocks are rebuilt
    for src, dst in ((FRAMES, 'frames_built.npz'), ('data/tethers.json', 'data/tethers_built.json')):
        if not a.dry_run and not os.path.exists(os.path.join(ROOT, dst)):
            shutil.copy(os.path.join(ROOT, src), os.path.join(ROOT, dst))
    for p in range(1, 4):
        print(f'  repair pass {p}')
        run('relax_frames.py', '--frames', FRAMES, '--out', FRAMES, '--problems', 'qc_problems.csv', dry=a.dry_run)
        run('tethers_layer3.py', '--frames', FRAMES, dry=a.dry_run)
        if run('check_standin.py', '--frames', FRAMES, dry=a.dry_run, check=False) == 0:
            break
    run('pairing_scorecard.py', '--frames', FRAMES, dry=a.dry_run, check=False)
    run('run_stage2.py', '--frames', FRAMES, '--outdir', 'blocks_tethered', '--size', a.size,
        '--workers', workers(a, 24), '--adopt', 'frames_built.npz', '--adopt-tethers', 'data/tethers_built.json',
        dry=a.dry_run)
    qc(a)


def final(a):
    qc(a, 'blocks_final')
    if a.dry_run or problems() > 0:
        if not a.dry_run:          # keep this list: the repeated check below overwrites qc_problems.csv
            shutil.copy(os.path.join(ROOT, 'qc_problems.csv'), os.path.join(ROOT, 'qc_problems_final.csv'))
        run('repair_tails.py', '--indir', 'blocks_final', '--outdir', 'blocks_final', '--from-qc', 'qc_problems_final.csv',
            '--size', a.size, dry=a.dry_run)
        qc(a, 'blocks_final')
    print('  final check written to qc_summary.txt' + ('' if a.dry_run else f' ({problems()} problems listed)'))


def options(ap):
    ap.add_argument('--size', type=int, default=50, help='nucleosomes per block (must match builder.py --size)')


if __name__ == '__main__':
    main('checker', __doc__, {'qc': qc, 'repair': repair, 'final': final}, options)
