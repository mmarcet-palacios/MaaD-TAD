"""ensemble_member.py -- build one member of the MaaD-TAD path ensemble (Folder + Placer).

Each member differs from the published model in three ways, all drawn from one seed:
  fold      a CSynth-like fold of the homology map (control_loci/PRAMEF_seqsim_like.csv; calibrated
            homology weight k = 0.007, homologous pairs only, as in Fig 8C) from a new random start
  spacing   every run of placeholder nucleosomes between two measured dyads is re-spaced at its own
            repeat length, drawn uniformly from 178-203 bp (the measured range); measured dyads are kept
  cohesin   each of the three tethers (T_I2, T_I5, T_I7) is present or absent with probability 0.5
Then Folder (refold_guide, make_naked, refold_tether) and Placer (stage1_frames, tethers_layer3,
relax_frames, tethers_layer3) run with the published settings.

    python ensemble_member.py --seed 7 [--out ensemble]     -> ensemble/m007/frames.npz, member.json
"""
import argparse, json, os, re, subprocess, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'control_loci'))
from csynth_like import fold

N_SEC = 325
LOCUS_LO, LOCUS_HI = 12771345, 13423345
NRL_LO, NRL_HI = 178, 203


def run(script, *args, cwd, log):
    t0 = time.time()
    cmd = [sys.executable, os.path.join(ROOT, script)] + [str(a) for a in args]
    with open(log, 'a') as f:
        f.write('\n> ' + ' '.join([script] + [str(a) for a in args]) + '\n'); f.flush()
        rc = subprocess.call(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    if rc:
        raise SystemExit(f'{script} failed (exit {rc}); see {log}')
    return round(time.time() - t0, 1)


def last_domain(log):
    """The --domain value the previous step printed for the next one (keeps the path at its own scale)."""
    return float(re.findall(r'--domain\s+([0-9.]+)', open(log).read().split('\n> ')[-1])[-1])


def member_fold(seed, out, k=0.007):
    S = np.loadtxt(os.path.join(ROOT, 'control_loci', 'PRAMEF_seqsim_like.csv'), delimiter=',', skiprows=1,
                   converters={0: lambda s: 0}, usecols=range(1, N_SEC + 1))[:N_SEC, :N_SEC]
    # 326 beads (section boundaries), as in the CSynth path: repeat the last section
    M = np.full((N_SEC + 1, N_SEC + 1), np.nan); M[:N_SEC, :N_SEC] = S
    M[N_SEC, :N_SEC] = S[-1]; M[:N_SEC, N_SEC] = S[:, -1]; M[N_SEC, N_SEC] = 100
    X, e = fold(M, kc=k, rep=0.3, seeds=1, rng0=10000 + seed, thresh=1.1)
    ids = 1 + 2000 * np.arange(N_SEC + 1)
    np.savetxt(out, np.c_[ids, X], fmt=['%d', '%.3f', '%.3f', '%.3f'], delimiter='\t')
    return float(e)


def member_bed(rng, src, out):
    head, rows = [], []
    for l in open(src):
        if l.startswith('#'):
            head.append(l)
        else:
            rows.append(l.rstrip('\n').split('\t'))
    meas = [r for r in rows if r[7] == 'real']
    anchors = [LOCUS_LO - 1 + 74 - 190] + [int(r[6]) for r in meas] + [LOCUS_HI + 1 - 74 + 190]
    # regions of the original placeholders, by position, to keep the region column
    ph = np.array([int(r[6]) for r in rows if r[7] == 'placeholder']); ph_reg = [r[8] for r in rows if r[7] == 'placeholder']
    new, runs = [], []
    mi = 0
    for a, b in zip(anchors[:-1], anchors[1:]):
        if mi > 0:
            new.append(meas[mi - 1])
        mi += 1
        span = b - a
        nrl = rng.uniform(NRL_LO, NRL_HI)
        k = int(round(span / nrl)) - 1
        if k < 1:
            continue
        step = span / (k + 1)
        if step < 163:             # keep the published floor
            k = int(span // 163) - 1; step = span / (k + 1)
        runs.append((a, b, nrl, k))
        for j in range(1, k + 1):
            c = int(round(a + j * step))
            if c - 73 < LOCUS_LO or c + 73 > LOCUS_HI:
                continue
            reg = ph_reg[int(np.argmin(np.abs(ph - c)))] if len(ph) else 'mapped'
            new.append(['chr1', str(c - 74), str(c + 73), 'gapfill', '.', '.', str(c), 'placeholder', reg])
    new.sort(key=lambda r: int(r[6]))
    for i, r in enumerate(new):
        if r[7] == 'placeholder':
            r[3] = f'gapfill_{i + 1}'
    nm = sum(r[7] == 'real' for r in new)
    with open(out, 'w') as f:
        f.write(head[0])
        f.write(f'# {len(new)} nucleosomes: {nm} measured, {len(new) - nm} placeholder (ensemble re-spacing {NRL_LO}-{NRL_HI} bp)\n')
        f.writelines(head[2:])
        for r in new:
            f.write('\t'.join(r) + '\n')
    d = np.array([int(r[6]) for r in new])
    return dict(n_nuc=len(new), n_measured=nm, min_nrl=int(np.diff(d).min()), median_nrl=float(np.median(np.diff(d))),
                run_nrls=[round(r[2], 1) for r in runs])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--out', default=os.path.join(ROOT, 'ensemble'))
    ap.add_argument('--p-tether', type=float, default=0.5)
    ap.add_argument('--k', type=float, default=0.007, help='homology weight of the fold (0.007 = calibrated on Micro-C, Fig 8C)')
    ap.add_argument('--no-relax', action='store_true', help='stop after stage 1 and the ring geometry')
    a = ap.parse_args()
    d = os.path.join(a.out, f'm{a.seed:03d}'); os.makedirs(d, exist_ok=True)
    log = os.path.join(d, 'log.txt'); open(log, 'w').close()
    rng = np.random.default_rng(a.seed)
    info = dict(seed=a.seed, k=a.k, times={})
    t0 = time.time()
    info['fold_energy'] = member_fold(a.seed, os.path.join(d, 'fold.xyz'), a.k); info['times']['fold'] = round(time.time() - t0, 1)
    info.update(member_bed(rng, os.path.join(ROOT, 'data', 'nucleosomes_GRCh38_resolved.bed'), os.path.join(d, 'nuc.bed')))
    R = json.load(open(os.path.join(ROOT, 'data', 'homology_restraints.json')))
    keep = [t for t in R['tethers'] if rng.random() < a.p_tether]
    R['tethers'] = keep; info['tethers'] = [t['id'] for t in keep]
    json.dump(R, open(os.path.join(d, 'restraints.json'), 'w'), indent=1)
    D = lambda n: os.path.join(ROOT, 'data', n)
    T = info['times']
    T['refold'] = run('refold_guide.py', '--bed', 'nuc.bed', '--guide', 'fold.xyz', '--gaps', D('gap_sections.txt'),
                      '--no-coil', '--coil-radius', 150, '--domain', 800, '--bend-radius', 1300, '--ref-min-sep', 6,
                      '--w-ref', 0.0005, '--out', 'guide_refolded.xyz', cwd=d, log=log)
    dom1 = last_domain(log)
    T['naked'] = run('make_naked.py', '--bed', 'nuc.bed', '--restraints', 'restraints.json', '--out', 'nuc_tethered.bed',
                     cwd=d, log=log)
    T['tether'] = run('refold_tether.py', '--bed', 'nuc_tethered.bed', '--guide', 'guide_refolded.xyz',
                      '--gaps', D('gap_sections.txt'), '--domain', dom1, '--restraints', 'restraints.json', '--w-ref', 0.02,
                      '--out', 'guide_tethered.xyz', cwd=d, log=log)
    info['domain_refold'], info['domain_tether'] = dom1, last_domain(log)
    T['frames'] = run('stage1_frames.py', '--bed', 'nuc_tethered.bed', '--guide', 'guide_tethered.xyz',
                      '--gaps', D('gap_sections.txt'), '--domain', last_domain(log), '--coil-radius', 0, '--smooth', 1,
                      '--calib', 0, '--look', 4, '--out', 'frames_stage1.npz', cwd=d, log=log)
    T['rings0'] = run('tethers_layer3.py', '--frames', 'frames_stage1.npz', '--restraints', 'restraints.json',
                      '--out', 'tethers.json', cwd=d, log=log)
    if a.no_relax:
        info['total_s'] = round(time.time() - t0, 1)
        json.dump(info, open(os.path.join(d, 'member.json'), 'w'), indent=1)
        print(json.dumps(info)); return
    T['relax'] = run('relax_frames.py', '--frames', 'frames_stage1.npz', '--out', 'frames.npz', '--tethers', 'tethers.json',
                     cwd=d, log=log)
    T['rings'] = run('tethers_layer3.py', '--frames', 'frames.npz', '--restraints', 'restraints.json',
                     '--out', 'tethers.json', cwd=d, log=log)
    info['total_s'] = round(time.time() - t0, 1)
    json.dump(info, open(os.path.join(d, 'member.json'), 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
