"""placer.py -- MaaD-TAD module 2, Placer: place the nucleosomes along the refolded path.

Steps (run in this order; settings as in S2 Table):
  frames   stage1_frames.py      one frame (origin + rotation) per nucleosome on a two-start zigzag that follows
                                 data/guide_tethered.xyz, anchored to genomic position, 4 steps tested ahead
                                 -> frames_tethered.npz
  relax    relax_frames.py       removes collisions (minimum 12 nm between nucleosomes, linker extension <= 0.93)
  rings    tethers_layer3.py     cohesin ring geometry and waypoints for the tether linkers -> data/tethers.json
  score    pairing_scorecard.py  prints how closely each homology block and tether is paired (report only)

Run from the MaaD-TAD folder:
    python placer.py
    python placer.py --steps relax rings
"""
from maadtad_steps import run, main

FRAMES = 'frames_tethered.npz'


def frames(a):
    run('stage1_frames.py', '--bed', 'data/nucleosomes_tethered.bed', '--guide', 'data/guide_tethered.xyz',
        '--gaps', 'data/gap_sections.txt', '--domain', 1011.5, '--coil-radius', 0, '--smooth', 1, '--calib', 0,
        '--look', 4, '--out', FRAMES, dry=a.dry_run)


def relax(a):
    run('relax_frames.py', '--frames', FRAMES, '--out', FRAMES, dry=a.dry_run)


def rings(a):
    run('tethers_layer3.py', '--frames', FRAMES, dry=a.dry_run)


def score(a):
    run('pairing_scorecard.py', '--frames', FRAMES, dry=a.dry_run, check=False)


if __name__ == '__main__':
    main('placer', __doc__, {'frames': frames, 'relax': relax, 'rings': rings, 'score': score})
