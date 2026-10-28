"""score_models.py -- score the full-atom test models against Micro-C, step by step through the pipeline.

For each locus it compares three structures with H1 Micro-C (and HFFc6):
  CSynth   <name>_H1.xyz                        the CSynth fold (input)
  Folder   <name>_H1/guide_folded.xyz           after refold_guide.py (fibre-length bonds, stiffness, excluded volume)
  model    <name>_H1/frames.npz                 the built model: nucleosome centres, averaged per 2 kb bin
                                                (the frames the blocks were built from; tail repair does not move cores)
Scores (same as compare_folds.py):
  r 2 kb   distance-corrected agreement: at each separation, z-score of log Micro-C vs z-score of -log distance
  r 10 kb  same on 10 kb blocks (model contacts as mean distance^-3)
  lines    SEQSIM line pairs (SEQSIM > 1.1, >= 3 bins apart): median distance relative to the median at the
           same separation, and % closer than that median
  vs CSynth  distance-corrected correlation of the structure with the CSynth fold itself
Also: domain size (largest extent, radius of gyration, nm).

Run in MaaD-TAD\\hic_test:
    python score_models.py            -> prints the table, writes score_models_results.txt
"""
import glob, os, warnings
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = 2000
LOCI = {'PWRN_chr15': 24_130_000, 'IGHC_chr14': 105_540_000}


def read_csv(p):
    rows = [l.rstrip('\n').split(',') for l in open(p)]
    return np.array([[float(v) for v in r[1:]] for r in rows[1:]])


def read_xyz(p):
    return np.array([[float(v) for v in l.split()[1:4]] for l in open(p) if len(l.split()) >= 4])


def model_beads(p, start, n):
    F = np.load(p)
    O, pos = F['origins'], F['dyad_pos']
    k = np.clip((pos - start) // BIN, 0, n - 1)
    X = np.full((n, 3), np.nan)
    for b in range(n):
        if (k == b).any():
            X[b] = O[k == b].mean(0)
    return X, O


def dist(X):
    return np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))


def zsep(M, mask, smin):
    n = len(M); Z = np.full_like(M, np.nan)
    for d in range(smin, n - 2):
        i = np.arange(n - d); j = i + d; ok = mask[i, j]
        if ok.sum() < 5: continue
        v = M[i, j][ok]; sd = v.std()
        if sd == 0: continue
        Z[i[ok], j[ok]] = (v - v.mean()) / sd
    return Z


def r_2kb(D, C, good, smin=2):
    n = len(D); mask = np.zeros((n, n), bool); mask[np.triu_indices(n, smin)] = True
    mask &= good[:, None] & good[None, :] & ~np.isnan(D)
    zC = zsep(np.log(C + 1e-3 * np.nanmedian(C[mask])), mask, smin); zD = zsep(-np.log(D), mask, smin)
    ok = ~np.isnan(zC) & ~np.isnan(zD)
    return np.corrcoef(zC[ok], zD[ok])[0, 1]


def coarse(M, g=None, k=5, isD=False):
    n = len(M) // k * k; A = M[:n, :n].astype(float).copy()
    if isD:
        A = A ** -3.0; np.fill_diagonal(A, np.nan)
    else:
        A[~g[:n], :] = np.nan; A[:, ~g[:n]] = np.nan
    return np.nanmean(A.reshape(n // k, k, n // k, k), axis=(1, 3))


def r_10kb(A, B):
    m = len(A); msk = np.triu(np.ones((m, m), bool), 1) & ~np.isnan(A) & ~np.isnan(B)
    a = zsep(np.log(A), msk, 1); b = zsep(np.log(B), msk, 1); ok = ~np.isnan(a) & ~np.isnan(b)
    return np.corrcoef(a[ok], b[ok])[0, 1]


def lines_score(D, lines):
    n = len(D); rel = []
    for i, j in zip(*np.where(lines)):
        d = j - i; k = np.arange(n - d); med = np.nanmedian(D[k, k + d])
        rel.append(D[i, j] / med)
    rel = np.array(rel)
    return np.nanmedian(rel), np.nanmean(rel < 1)


def fold_r(D1, D2):
    n = len(D1); m = np.triu(np.ones((n, n), bool), 2) & ~np.isnan(D1) & ~np.isnan(D2)
    a = zsep(-np.log(D1), m, 2); b = zsep(-np.log(D2), m, 2); ok = ~np.isnan(a) & ~np.isnan(b)
    return np.corrcoef(a[ok], b[ok])[0, 1]


def size_nm(X, unit):
    X = X[~np.isnan(X[:, 0])]
    rg = np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean())
    return dist(X).max() * unit, rg * unit


def main():
    os.chdir(HERE); np.seterr(divide="ignore", invalid="ignore"); warnings.simplefilter("ignore", RuntimeWarning)
    out = []
    say = lambda s: (print(s), out.append(s))
    for name, start in LOCI.items():
        builds = sorted(d for d in glob.glob(f'{name}_H1*') if os.path.exists(os.path.join(d, 'frames.npz')))
        if not builds:
            say(f'== {name}: no built model yet, skipped'); continue
        z = np.load(f'{name}_microc.npz'); sim = read_csv(f'{name}_SEQSIM.csv'); n = len(sim)
        H1 = z['H1_bal_2000'].astype(float); HF = z['HFFc6_bal_2000'].astype(float)
        gH1 = ~np.isnan(z['H1_w_2000']); gHF = ~np.isnan(z['HFFc6_w_2000'])
        lines = np.triu(sim > 1.1, 3)
        Xc = read_xyz(f'{name}_H1.xyz')
        S = {'CSynth': (Xc, None)}
        for W in builds:
            tag = W[len(name) + 4:].lstrip('_') or 'cut all'
            Xm, O = model_beads(os.path.join(W, 'frames.npz'), start, n)
            S[f'Folder {tag}'] = (read_xyz(os.path.join(W, 'guide_folded.xyz')), 0.1)
            S[f'model {tag}'] = (Xm, 0.1)
        Dc = dist(Xc)
        cH = coarse(H1, gH1)
        say(f'== {name}: {n} bins of 2 kb, {len(O)} nucleosomes, {int(lines.sum())} line pairs')
        say(f'   {"":15s} {"r 2kb H1":>9s} {"r 2kb HFF":>10s} {"r 10kb H1":>10s} {"lines":>6s} {"closer":>7s} '
            f'{"vs CSynth":>10s} {"extent nm":>10s} {"Rg nm":>7s}')
        for lab, (X, unit) in S.items():
            D = dist(X)
            l_rel, l_cl = lines_score(D, lines)
            ext, rg = size_nm(X, unit) if unit else (np.nan, np.nan)
            say(f'   {lab:15s} {r_2kb(D, H1, gH1):+9.2f} {r_2kb(D, HF, gHF):+10.2f} '
                f'{r_10kb(coarse(D, isD=True), cH):+10.2f} {l_rel:6.2f} {100 * l_cl:6.0f}% '
                f'{(fold_r(D, Dc) if lab != "CSynth" else 1.0):+10.2f} {ext:10.0f} {rg:7.0f}')
        # nucleosome-level view: how far apart are the parallel lines in the built model(s)?
        i, j = np.where(lines)
        for lab, (X, unit) in S.items():
            if not lab.startswith('model'): continue
            dn = dist(X)
            say(f'   {lab}: line pairs median {np.nanmedian(dn[i, j]) / 10:.0f} nm apart '
            f'(all pairs at the same separations: {np.nanmedian([np.nanmedian(np.diagonal(dn, d)) for d in (j - i)]) / 10:.0f} nm)')
    open('score_models_results.txt', 'w').write('\n'.join(out) + '\n')
    print('wrote score_models_results.txt')


if __name__ == '__main__':
    main()
