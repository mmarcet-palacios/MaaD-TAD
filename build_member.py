"""build_member.py -- take one ensemble member (from ensemble_member.py --no-relax) to atoms.

Runs, inside the member folder:
  relax    Placer: relax_frames (clears colliding nucleosomes), then tethers_layer3
  reach    relax_frames --only-reach: moves only the nucleosomes around linkers that the relaxation left stretched
           past their contour length, until every linker is within 0.95 of it; then tethers_layer3.
           (A linker longer than its DNA cannot be closed: in member m030 the 24 linkers that Builder could not
           close were exactly the 24 stretched past their contour, and every other linker closed.)
  build    Builder: run_stage2 (all-atom blocks of 50 nucleosomes). Blocks whose frames changed since they were
           built are rebuilt; the others are kept.
  qc       Checker: stage4_qc. With --repair, up to 3 repair passes (relax_frames --problems, reach, tethers_layer3,
           check_standin), a rebuild of the changed blocks, and qc again.
  seq      optional (--seq): thread_sequence on the GRCh38 locus, repair_tails, final qc
           (skipped by default: it roughly doubles the disk per member)

    python build_member.py ensemble/m030 --workers 24                     # relax, reach, build, qc
    python build_member.py ensemble/m030 --workers 24 --steps reach build qc   # after an earlier relax and build
    python build_member.py ensemble/m030 --steps qc --repair              # qc with the repair loop
    python build_member.py ensemble/m030 --dry-run                        # print the commands only

Resumable: run_stage2 skips finished blocks, and relax is skipped if frames.npz exists.
All paths are passed absolute, because run_stage2 runs stage2_block from the repository folder.
Writes build.json (times, qc summary) in the member folder.
"""
import argparse, csv, glob, json, os, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
STEPS = ['relax', 'reach', 'build', 'qc', 'seq']
REACH = ['--only-reach', '--reach-weight', 1, '--reach-trigger', 0.95, '--win', 3, '--rounds', 3, '--maxiter', 500]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('member', help='member folder, e.g. ensemble/m030')
    ap.add_argument('--workers', type=int, default=24)
    ap.add_argument('--size', type=int, default=50)
    ap.add_argument('--steps', nargs='*', choices=STEPS, help='default: relax reach build qc (plus seq with --seq)')
    ap.add_argument('--repair', action='store_true', help='qc step: run the repair loop if problems are found')
    ap.add_argument('--seq', action='store_true', help='also thread the sequence and repair tails')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    d = os.path.abspath(a.member)
    if not os.path.exists(os.path.join(d, 'frames_stage1.npz')):
        sys.exit(f'{d} has no frames_stage1.npz: run ensemble_member.py first')
    steps = a.steps or (['relax', 'reach', 'build', 'qc'] + (['seq'] if a.seq else []))
    P = lambda n: os.path.join(d, n)
    log = P('build_log.txt')
    info = json.load(open(P('build.json'))) if os.path.exists(P('build.json')) else dict(times={})

    def run(script, *args, check=True):
        cmd = [sys.executable, os.path.join(ROOT, script)] + [str(x) for x in args]
        print('>', ' '.join([script] + [str(x) for x in args]), flush=True)
        if a.dry_run:
            return 0
        t0 = time.time()
        with open(log, 'a') as f:
            f.write('\n> ' + ' '.join(cmd) + '\n'); f.flush()
            rc = subprocess.call(cmd, cwd=d, stdout=f, stderr=subprocess.STDOUT)
        info['times'][script] = round(info['times'].get(script, 0) + time.time() - t0, 1)
        if rc and check:
            sys.exit(f'{script} failed (exit {rc}); see {log}')
        return rc

    def rings():
        run('tethers_layer3.py', '--frames', P('frames.npz'), '--restraints', P('restraints.json'), '--out', P('tethers.json'))

    def reach():
        run('relax_frames.py', '--frames', P('frames.npz'), '--out', P('frames.npz'), '--tethers', P('tethers.json'), *REACH)

    def open_linkers(blocks):
        """Linkers Builder could not close, in blocks that have no finished .cif (block_XX_FAILED.cif only)."""
        out = []
        for f in sorted(glob.glob(os.path.join(P(blocks), 'block_??_FAILED.cif'))):
            q = f.replace('_FAILED.cif', '_qc.csv')
            if os.path.exists(q) and not os.path.exists(f.replace('_FAILED.cif', '.cif')):
                out += [int(r['linker']) for r in csv.DictReader(open(q)) if r.get('ok') != 'True']
        return out

    def qc(blocks):
        run('stage4_qc.py', '--frames', P('frames.npz'), '--blocks', P(blocks), '--size', a.size,
            '--workers', min(a.workers, 16))
        if a.dry_run:
            return 0
        rows = list(csv.DictReader(open(P('qc_problems.csv'))))
        # blocks with an open linker are left out of the model, so stage4_qc cannot see them; list their open
        # linkers as breaks so that the repair pass relaxes the frames there
        failed = open_linkers(blocks)
        if failed:
            print(f'  open linkers in failed blocks: {failed}')
            with open(P('qc_problems.csv'), 'a', newline='') as fh:
                w = csv.writer(fh)
                for L in failed:
                    w.writerow(['break', '', '', L, '', '', ''])
        info['failed_linkers'] = failed
        return len(rows) + len(failed)

    def stage2(*extra):
        run('run_stage2.py', '--frames', P('frames.npz'), '--outdir', P('blocks_tethered'), '--size', a.size,
            '--workers', a.workers, '--tethers', P('tethers.json'), *extra)
        if a.dry_run:
            return
        # a block that failed earlier and has now been built: drop its old *_FAILED.cif
        for f in glob.glob(os.path.join(P('blocks_tethered'), 'block_??_FAILED.cif')):
            if os.path.exists(f.replace('_FAILED.cif', '.cif')):
                os.remove(f)
        info['blocks_built'] = len(glob.glob(os.path.join(P('blocks_tethered'), 'block_??.cif')))
        info['open_linkers'] = open_linkers('blocks_tethered')
        print(f"  blocks built {info['blocks_built']}; open linkers {info['open_linkers']}")

    if 'relax' in steps and not os.path.exists(P('frames.npz')):
        run('relax_frames.py', '--frames', P('frames_stage1.npz'), '--out', P('frames.npz'), '--tethers', P('tethers.json'),
            '--reach-weight', 1)
        rings()
    if 'reach' in steps:
        reach(); rings()
    if 'build' in steps:
        stage2()
    if 'qc' in steps:
        n = qc('blocks_tethered'); info['qc_problems_first'] = info.get('qc_problems_first', n)
        if n and a.repair:
            for src, dst in (('frames.npz', 'frames_built.npz'), ('tethers.json', 'tethers_built.json')):
                if not a.dry_run and not os.path.exists(P(dst)):
                    shutil.copy(P(src), P(dst))
            for p in range(1, 4):
                print(f'  repair pass {p}')
                run('relax_frames.py', '--frames', P('frames.npz'), '--out', P('frames.npz'),
                    '--problems', P('qc_problems.csv'), '--tethers', P('tethers.json'), '--reach-weight', 1)
                reach(); rings()
                if run('check_standin.py', '--frames', P('frames.npz'), '--tethers', P('tethers.json'), check=False) == 0:
                    break
            info['repair_passes'] = info.get('repair_passes', 0) + p
            stage2('--adopt', P('frames_built.npz'), '--adopt-tethers', P('tethers_built.json'))
            n = qc('blocks_tethered')
        info['qc_problems'] = n
        if not a.dry_run:
            info['qc_summary'] = open(P('qc_summary.txt')).read().splitlines()
    if 'seq' in steps:
        run('thread_sequence.py', '--indir', P('blocks_tethered'), '--outdir', P('blocks_seq'),
            '--fasta', os.path.join(ROOT, 'locus_GRCh38.fa'), '--workers', min(a.workers, 6))
        run('repair_tails.py', '--indir', P('blocks_seq'), '--outdir', P('blocks_final'), '--size', a.size,
            '--workers', min(a.workers, 6))
        info['qc_problems_final'] = qc('blocks_final')
        if not a.dry_run:
            info['qc_summary_final'] = open(P('qc_summary.txt')).read().splitlines()
    if not a.dry_run:
        json.dump(info, open(P('build.json'), 'w'), indent=1)
        print(json.dumps({k: v for k, v in info.items() if not k.startswith('qc_summary')}))


if __name__ == '__main__':
    main()
