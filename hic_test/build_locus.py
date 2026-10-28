"""build_locus.py -- full-atom MaaD-TAD build of a test locus from its CSynth fold (Micro-C only: no homology, no tethers).

Same pipeline and settings as the PRAMEF model, minus the pairing layer:
  bed     placeholder nucleosome ladder (no measured positions here): repeat 190-191 bp, chosen so the
          chain spans the window exactly (first core starts at the first base, last core ends at the last)
  fold    Folder: refold_guide.py --no-coil on <name>_<fold>.xyz (bend radius 130 nm, excluded volume
          52 nm, w-ref 0.02, as S2 Table). CSynth units are scaled so the mean bead spacing is 0.62 x the
          fibre length each 2 kb bin's nucleosomes supply: the same starting density as PRAMEF (800 nm start).
  frames  Placer: stage1_frames.py (--coil-radius 0 --smooth 1 --calib 0 --look 4, as PRAMEF)
  relax   relax_frames.py (no tethers)
  build   Builder: run_stage2.py (blocks of 50 nucleosomes)
  qc      Checker: stage4_qc.py on blocks/
  seq     thread_locus.py: GRCh38 bases from <name>.fa
  tails   repair_tails.py: blocks_seq/ -> blocks_final/
  final   stage4_qc.py on blocks_final/ (writes into qc_final/)
  repair  (only if qc lists problems) relax with qc_problems.csv, rebuild the changed blocks, qc again;
          then rerun seq, tails, final

Everything for one locus goes in hic_test\\<name>_<fold>\\ (e.g. PWRN_chr15_H1\\); PRAMEF files are not touched.
The PRAMEF scripts are used unchanged from MaaD-TAD; the locus start replaces PRAMEF's
ROI_START inside this process only.

Run in MaaD-TAD\\hic_test, chromatin env:
    python build_locus.py PWRN_chr15 IGHC_chr14              # both, one after the other (all steps)
    python build_locus.py IGHC_chr14 --steps seq tails final # some steps only
    python build_locus.py PWRN_chr15 --steps repair
Options: --fold H1|hybrid|SEQSIM|HFFc6 (default H1), --workers 24, --domain <nm> (fix the start scale),
         --ref-cut <nm>: Folder keeps only CSynth distances shorter than this (PRAMEF: 300); default 0 = keep all.
         A cut gets its own folder, e.g. --ref-cut 300 -> PWRN_chr15_H1_cut300/
"""
import argparse, glob, json, os, shutil, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.dirname(HERE)
LOCI = {'PWRN_chr15': ('chr15', 24_130_000, 24_440_000),
        'IGHC_chr14': ('chr14', 105_540_000, 105_750_000)}
BIN, CORE, NRL = 2000, 147, 190
ALL = ['bed', 'fold', 'frames', 'relax', 'build', 'qc', 'seq', 'tails', 'final']
FOLD_ARGS = ['--coil-radius', '150', '--bend-radius', '1300', '--w-ref', '0.02', '--no-coil']
START_RATIO = 0.62          # PRAMEF: mean bead spacing at --domain 800 / fibre per bin
FRAME_ARGS = ['--coil-radius', '0', '--smooth', '1', '--calib', '0', '--look', '4']


def read_xyz(path):
    return np.array([[float(v) for v in l.split()[1:4]] for l in open(path) if len(l.split()) >= 4])


def write_bed(path, chrom, start, end):
    """Placeholder ladder spanning start+1..end (1-based): dyad i at start+74+sum(nrl[:i])."""
    L = end - start - CORE
    n = L // NRL                                  # number of linkers
    extra = L - n * NRL                           # spread the remainder as +1 bp, evenly
    nrl = NRL + (np.floor((np.arange(n) + 1) * extra / n) - np.floor(np.arange(n) * extra / n)).astype(int)
    c = start + 74 + np.concatenate([[0], np.cumsum(nrl)])
    assert c[-1] + 73 == end and c[0] - 73 == start + 1
    with open(path, 'w') as f:
        f.write(f'# GRCh38 {chrom}:{start + 1}-{end}\n')
        f.write(f'# {len(c)} nucleosomes: 0 measured, {len(c)} placeholder\n')
        f.write(f'# placeholder ladder, NRL {nrl.min()}-{nrl.max()} bp (mean {nrl.mean():.2f})\n')
        f.write('# chrom\tstart\tend\tname\tscore\tstrand\tcenter\tsource\tregion\n')
        for i, x in enumerate(c):
            f.write(f'{chrom}\t{x - 74}\t{x + 73}\tplaceholder_{i + 1}\t.\t.\t{x}\tplaceholder\tmapped\n')
    return len(c)


def patched(start):
    """PRAMEF modules with ROI_START = this locus (bead k of the guide = bin start + k*2 kb)."""
    sys.path.insert(0, MODEL)
    import stage1_frames as s1
    s1.ROI_START, s1.FIRST_SECTION = start, 1
    import refold_guide as rg
    rg.ROI_START, rg.FIRST_SECTION = start, 1
    return s1, rg


def call(argv, s1_or_rg):
    old = sys.argv; sys.argv = argv
    try:
        s1_or_rg.main()
    finally:
        sys.argv = old


def run(script, args, cwd):
    cmd = [sys.executable, script if os.path.isabs(script) else os.path.join(MODEL, script)] + args
    print('  >', ' '.join(os.path.basename(c) if os.path.isabs(c) and c.endswith('.py') else c for c in cmd[1:]), flush=True)
    rc = subprocess.call(cmd, cwd=cwd)
    if rc:
        sys.exit(f'{script} failed (exit {rc})')


def problems(path):
    return max(0, sum(1 for _ in open(path)) - 1) if os.path.exists(path) else -1


def build(name, fold, steps, workers, domain, ref_cut=0):
    chrom, start, end = LOCI[name]
    W = os.path.join(HERE, f'{name}_{fold}' + (f'_cut{ref_cut:g}' if ref_cut else ''))
    os.makedirs(W, exist_ok=True)
    P = lambda *p: os.path.join(W, *p)
    bed, gaps, guide, frames = P('nucleosomes.bed'), P('gaps.txt'), P('guide_folded.xyz'), P('frames.npz')
    info_p = P('build_info.json')
    info = json.load(open(info_p)) if os.path.exists(info_p) else {}
    log = open(P('build_log.txt'), 'a')
    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()
    say(f'\n######## {name}  {chrom}:{start + 1:,}-{end:,}  fold {fold}  -> {os.path.basename(W)}  ({time.strftime("%Y-%m-%d %H:%M")})')
    s1, rg = patched(start)
    for step in steps:
        t0 = time.time()
        say(f'=== {step}')
        if step == 'bed':
            n = write_bed(bed, chrom, start, end)
            open(gaps, 'w').close()
            info.update(nucleosomes=n); say(f'  {n} placeholder nucleosomes -> nucleosomes.bed; gaps.txt empty')
        elif step == 'fold':
            src = os.path.join(HERE, f'{name}_{fold}.xyz')
            x = read_xyz(src)
            assert len(x) == (end - start) // BIN, f'{src}: {len(x)} beads, expected {(end - start) // BIN}'
            need = rg.demand(bed, len(x) - 1)                       # fibre axis per bin, A
            bond = np.linalg.norm(np.diff(x, axis=0), axis=1)
            dmax = np.linalg.norm(x[:, None] - x[None], axis=-1).max()
            d0 = domain or dmax * START_RATIO * need.mean() / bond.mean() / 10.0  # nm
            say(f'  {len(x)} beads; fibre per 2 kb bin {need.mean() / 10:.1f} nm; start domain {d0:.1f} nm')
            call(['refold_guide.py', '--bed', bed, '--guide', src, '--gaps', gaps, '--domain', f'{d0:.2f}',
                  '--out', guide, '--ref-cut', f'{(ref_cut or 1e7) * 10:g}'] + FOLD_ARGS, rg)
            info.update(ref_cut_nm=ref_cut or 'all')
            g = read_xyz(guide)
            info.update(start_domain_nm=round(d0, 1),
                        domain_nm=round(np.linalg.norm(g[:, None] - g[None], axis=-1).max() / 10, 1))
            say(f'  folded guide: largest extent {info["domain_nm"]} nm')
        elif step == 'frames':
            call(['stage1_frames.py', '--bed', bed, '--guide', guide, '--gaps', gaps,
                  '--domain', str(info['domain_nm']), '--out', frames] + FRAME_ARGS, s1)
        elif step == 'relax':
            run('relax_frames.py', ['--frames', frames, '--out', frames, '--tethers', 'none'], W)
        elif step == 'build':
            run('run_stage2.py', ['--frames', frames, '--outdir', P('blocks'), '--workers', str(workers),
                                  '--tethers', 'none'], W)
        elif step == 'qc':
            run('stage4_qc.py', ['--frames', frames, '--blocks', P('blocks'), '--workers', '16'], W)
            say(f'  qc_problems.csv: {problems(P("qc_problems.csv"))} problems')
        elif step == 'seq':
            run(os.path.join(HERE, 'thread_locus.py'), ['--indir', P('blocks'), '--outdir', P('blocks_seq'),
                '--fasta', os.path.join(HERE, f'{name}.fa'), '--start', str(start), '--end', str(end),
                '--workers', '12'], W)
        elif step == 'tails':
            run('repair_tails.py', ['--indir', P('blocks_seq'), '--outdir', P('blocks_final'), '--workers', '12'], W)
        elif step == 'final':
            os.makedirs(P('qc_final'), exist_ok=True)
            # chain A's seam phosphate (first residue of each block) is stored in blocks/block_XX_seam.json, not
            # in the .cif; threading and tail repair leave the backbone untouched, so the same seams apply
            for f in glob.glob(P('blocks', 'block_??_seam.json')):
                shutil.copy2(f, P('blocks_final'))
            run('stage4_qc.py', ['--frames', frames, '--blocks', P('blocks_final'), '--workers', '16'], P('qc_final'))
            k = problems(P('qc_final', 'qc_problems.csv'))
            fails = [l.strip() for l in open(P('qc_final', 'qc_summary.txt')) if 'FAIL' in l]
            say(f'  final model: {k} steric problems, {len(fails)} failed checks' + ('' if k == 0 and not fails else
                f'  -> python build_locus.py {name} --fold {fold} --steps repair'))
            for l in fails:
                say('    ' + l)
        elif step == 'repair':
            if problems(P('qc_problems.csv')) <= 0:
                say('  qc_problems.csv lists nothing to repair'); continue
            run('relax_frames.py', ['--frames', frames, '--out', frames, '--problems', P('qc_problems.csv'),
                                    '--tethers', 'none'], W)
            run('run_stage2.py', ['--frames', frames, '--outdir', P('blocks'), '--workers', str(workers),
                                  '--tethers', 'none'], W)          # rebuilds only blocks whose frames changed
            run('stage4_qc.py', ['--frames', frames, '--blocks', P('blocks'), '--workers', '16'], W)
            say(f'  after repair: {problems(P("qc_problems.csv"))} problems '
                f'(then: --steps seq tails final; repeat repair if still > 0)')
        info[f'{step}_min'] = round((time.time() - t0) / 60, 1)
        json.dump(info, open(info_p, 'w'), indent=1)
        say(f'  {step} done in {info[f"{step}_min"]} min')
    log.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('loci', nargs='+', choices=list(LOCI))
    ap.add_argument('--fold', default='H1', choices=['H1', 'hybrid', 'SEQSIM', 'HFFc6'])
    ap.add_argument('--steps', nargs='+', default=ALL, choices=ALL + ['repair'])
    ap.add_argument('--workers', type=int, default=24)
    ap.add_argument('--domain', type=float, help='start scale for the Folder, nm (default: from fibre length)')
    ap.add_argument('--ref-cut', type=float, default=0, help='Folder keeps CSynth distances below this, nm (0 = all)')
    a = ap.parse_args()
    for name in a.loci:
        build(name, a.fold, a.steps, a.workers, a.domain, a.ref_cut)


if __name__ == '__main__':
    main()
