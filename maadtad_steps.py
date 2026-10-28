"""maadtad_steps.py -- shared helper for the four MaaD-TAD module drivers (folder.py, placer.py, builder.py, checker.py).

Each driver is a list of named steps; each step runs one or more of the pipeline scripts in this folder with the
settings used for the PRAMEF/HNRNPCL model (S2 Table). The scripts can also be run on their own (see README.md).
"""
import argparse, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))


def run(script, *args, dry=False, check=True):
    """Run one pipeline script from the repository root; returns its exit code."""
    cmd = [sys.executable, os.path.join(ROOT, script)] + [str(a) for a in args]
    print('  > python ' + ' '.join([script] + [str(a) for a in args]), flush=True)
    if dry:
        return 0
    rc = subprocess.call(cmd, cwd=ROOT)
    if check and rc:
        sys.exit(f'{script} failed (exit {rc})')
    return rc


def problems(path='qc_problems.csv'):
    """Number of problems listed by stage4_qc.py (-1 if the file is missing)."""
    p = os.path.join(ROOT, path)
    return max(0, sum(1 for _ in open(p)) - 1) if os.path.exists(p) else -1


def workers(a, paper):
    """Worker count: --workers if given, else the value used for the paper (24-core workstation), capped at the CPU count."""
    return a.workers or max(1, min(paper, os.cpu_count() or 1))


def main(name, doc, steps, extra=None):
    """steps: ordered dict {step name: function(args)}."""
    ap = argparse.ArgumentParser(prog=f'{name}.py', description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--steps', nargs='+', choices=list(steps), default=list(steps), help='steps to run (default: all, in order)')
    ap.add_argument('--dry-run', action='store_true', help='print the commands without running them')
    ap.add_argument('--workers', type=int, help='parallel processes (default: the value used for the paper, capped at the CPU count)')
    if extra:
        extra(ap)
    a = ap.parse_args()
    t0 = time.time()
    for s in steps:
        if s in a.steps:
            print(f'=== {name}: {s}', flush=True)
            steps[s](a)
    print(f'=== {name}: done ({time.time() - t0:.0f} s)')
