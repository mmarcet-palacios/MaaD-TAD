"""homology_test_loci.py -- the Fig 8 contact enrichment test, applied to the two Micro-C-only test loci (IGHC, PWRN).

Same test as compare_hic.py (control_loci, Fig 8B):
  O/E = balanced contacts divided by the mean of their diagonal (2 kb bins);
  homologous pairs: SEQSIM > 1.1, separated by >= 3 bins, both bins covered (finite, non-zero O/E);
  one control per homologous pair, drawn at random (seed 0) from pairs at the same separation with SEQSIM <= 0.6;
  one-sided Mann-Whitney U on the O/E values; enrichment = median of the per-pair ratios.
SEQSIM here is the real SEQSIM map of each window (<name>_SEQSIM.csv), not the SEQSIM-like map of the control loci.
Run on H1 (the map the models were built from), HFFc6 (a second cell type), and the full-atom model
(contacts = distance^-3 between mean nucleosome centres of each 2 kb bin, as in Fig 10B), which was built without homology.
Seed check: the median ratio over 200 seeds is printed beside the seed-0 value.

Run in MaaD-TAD\\hic_test (needs numpy, scipy):
    python homology_test_loci.py      -> homology_test_loci_results.txt
"""
import os
import numpy as np
from scipy.stats import mannwhitneyu
from fig_models_maps import HERE, BIN, LOCI, model_contacts
from fig_seqsim_maps import read_csynth


def oe(M):
    O = np.full_like(M, np.nan); n = len(M)
    for d in range(n):
        v = np.diagonal(M, d); e = np.nanmean(v) if np.isfinite(v).any() else np.nan
        idx = np.arange(n - d); O[idx, idx + d] = v / e if e else np.nan; O[idx + d, idx] = O[idx, idx + d]
    return O


def test(M, S, seed=0):
    n = len(M); Mo = oe(M)
    iu = np.triu_indices(n, 3); hv, sv = Mo[iu], S[iu]; sep = iu[1] - iu[0]
    ok = np.isfinite(hv) & (hv > 0)
    hom = ok & (sv > 1.1); rng = np.random.default_rng(seed); h, c = [], []
    for k in np.where(hom)[0]:
        cand = np.where(ok & (sv <= 0.6) & (sep == sep[k]))[0]
        if len(cand):
            j = rng.choice(cand); h.append(hv[k]); c.append(hv[j])
    h, c = np.array(h), np.array(c)
    return dict(pairs=len(h), med_h=np.median(h), med_c=np.median(c), ratio=np.median(h / c),
                p=mannwhitneyu(h, c, alternative='greater').pvalue)


if __name__ == '__main__':
    lines = ['Fig 8 enrichment test on the Micro-C-only test loci (2 kb bins; SEQSIM > 1.1 vs <= 0.6, same separation, >= 3 bins)',
             f'{"locus":6} {"map":8} {"pairs":>5} {"O/E hom":>8} {"O/E ctrl":>8} {"ratio":>6} {"p (1-sided)":>12} {"ratio, 200 seeds (median, 2.5-97.5%)":>40}']
    for tag, (name, chrom, start, end, _) in LOCI.items():
        n = (end - start) // BIN
        S = read_csynth(os.path.join(HERE, f'{name}_SEQSIM.csv'))
        z = np.load(os.path.join(HERE, f'{name}_microc.npz'))
        maps = {}
        for cell in ('H1', 'HFFc6'):
            M = z[f'{cell}_bal_2000'].astype(float); w = z[f'{cell}_w_2000']
            M[np.isnan(w), :] = np.nan; M[:, np.isnan(w)] = np.nan
            maps[cell] = M
        C = model_contacts(os.path.join(HERE, f'{name}_H1', 'frames.npz'), start, n)
        maps['model'] = np.where(np.isfinite(C), C, np.nan)
        for lab, M in maps.items():
            r = test(M, S)
            rs = np.array([test(M, S, s)['ratio'] for s in range(200)])
            lines.append(f'{tag:6} {lab:8} {r["pairs"]:5d} {r["med_h"]:8.2f} {r["med_c"]:8.2f} {r["ratio"]:6.2f} {r["p"]:12.2e} '
                         f'{np.median(rs):18.2f} ({np.percentile(rs, 2.5):.2f}-{np.percentile(rs, 97.5):.2f})')
    txt = '\n'.join(lines); print(txt)
    open(os.path.join(HERE, 'homology_test_loci_results.txt'), 'w').write(txt + '\n')
