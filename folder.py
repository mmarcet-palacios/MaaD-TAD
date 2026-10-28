"""folder.py -- MaaD-TAD module 1, Folder: refold the CSynth path so that it can hold a chromatin fibre,
then apply the homology pairing restraints and the cohesin tether terms.

Steps (run in this order; settings as in S2 Table):
  refold   refold_guide.py   CSynth path (data/gap-filled.xyz, 326 beads of 2 kb) -> data/guide_refolded.xyz
                             no coil, excluded volume 52 nm, minimum bend radius 130 nm, start at 800 nm
  naked    make_naked.py     opens at least 400 bp of naked DNA at each tethered site
                             data/nucleosomes_GRCh38_resolved.bed -> data/nucleosomes_tethered.bed (17 nucleosomes removed)
  tether   refold_tether.py  pairing (arm spacing, direction, register) and tether terms on the refolded path,
                             fidelity weight 0.02 -> data/guide_tethered.xyz (about 7 min)

Run from the MaaD-TAD folder, in the conda environment of environment.yml:
    python folder.py                    # all three steps
    python folder.py --steps tether     # one step
    python folder.py --dry-run          # print the commands only
The outputs used for the paper are supplied in data/, so Placer can be run without rerunning Folder.
"""
from maadtad_steps import run, main


def refold(a):
    run('refold_guide.py', '--bed', 'data/nucleosomes_GRCh38_resolved.bed', '--guide', 'data/gap-filled.xyz',
        '--gaps', 'data/gap_sections.txt', '--no-coil', '--coil-radius', 150, '--domain', 800, '--bend-radius', 1300,
        '--ref-min-sep', 6, '--w-ref', 0.0005, '--out', 'data/guide_refolded.xyz', dry=a.dry_run)


def naked(a):
    run('make_naked.py', dry=a.dry_run)


def tether(a):
    run('refold_tether.py', '--bed', 'data/nucleosomes_tethered.bed', '--guide', 'data/guide_refolded.xyz',
        '--w-ref', 0.02, '--out', 'data/guide_tethered.xyz', dry=a.dry_run)


if __name__ == '__main__':
    main('folder', __doc__, {'refold': refold, 'naked': naked, 'tether': tether})
