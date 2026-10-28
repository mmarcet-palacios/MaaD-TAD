"""to_contacts.py -- write the hic_test maps as CSynth .contacts files (3 columns: pos1, pos2, value; tab-separated).

Run in MaaD-TAD\\hic_test after extract_hic_test.py (needs numpy):
    python to_contacts.py
For each locus (PWRN_chr15, IGHC_chr14) it writes:
  <name>_H1.contacts       H1 Micro-C, balanced 2 kb, scaled as in <name>_H1_microc.csv; discarded bins left out
  <name>_HFFc6.contacts    same for HFFc6
  <name>_SEQSIM.contacts   SEQSIM map (max of direct and inverted), background 0.4
  <name>_hybrid.contacts   H1 Micro-C, with the discarded bins filled from SEQSIM. SEQSIM values are converted to
                           Micro-C values by rank, separately at each genomic separation (so the fill keeps the
                           distance decay of Micro-C), using the Micro-C values at that separation as the target.
Positions are genomic bin starts (GRCh38, multiples of 2 kb), one chromosome per file; only pairs i < j are written
(no diagonal), and pairs with value 0 are omitted. CSynth places a bead in every 2 kb bin of the range.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOCI = {'PWRN_chr15': 24_130_000, 'IGHC_chr14': 105_540_000}
BIN = 2000


def read_csv(path):
    rows = [ln.rstrip('\n').split(',') for ln in open(path)]
    return np.array([[float(v) for v in r[1:]] for r in rows[1:]])


def write_contacts(path, M, start, keep=None):
    n = len(M)
    k = 0
    with open(path, 'w') as f:
        for i in range(n):
            for j in range(i + 1, n):
                v = M[i, j]
                if v > 0 and (keep is None or (keep[i] and keep[j])):
                    f.write(f'{start + i * BIN}\t{start + j * BIN}\t{v:.5f}\n'); k += 1
    return k


def hybrid(mc, sim, disc):
    """Fill rows/cols of discarded bins in mc with SEQSIM, rank-matched to Micro-C at each separation."""
    n = len(mc)
    H = mc.copy()
    for d in range(1, n):
        i = np.arange(n - d); j = i + d
        good = ~(disc[i] | disc[j])
        fill = ~good
        if not fill.any():
            continue
        ref = np.sort(mc[i[good], j[good]]) if good.any() else None
        if ref is None or len(ref) == 0:
            continue
        s_all = sim[i, j]
        # rank of each SEQSIM value among all SEQSIM values at this separation -> same quantile of Micro-C
        q = (np.argsort(np.argsort(s_all)) + 0.5) / len(s_all)
        vals = np.quantile(ref, q[fill])
        H[i[fill], j[fill]] = vals; H[j[fill], i[fill]] = vals
    return H


def main():
    os.chdir(HERE)
    for name, start in LOCI.items():
        z = np.load(f'{name}_microc.npz')
        sim = read_csv(f'{name}_SEQSIM.csv')
        for cell in ('H1', 'HFFc6'):
            mc = read_csv(f'{name}_{cell}_microc.csv')
            disc = np.isnan(z[f'{cell}_w_2000'])
            k = write_contacts(f'{name}_{cell}.contacts', mc, start, keep=~disc)
            print(f'{name}_{cell}.contacts: {k} pairs ({disc.sum()} discarded bins left out)')
            if cell == 'H1':
                H = hybrid(mc, sim, disc)
                k = write_contacts(f'{name}_hybrid.contacts', H, start)
                print(f'{name}_hybrid.contacts: {k} pairs ({disc.sum()} bins filled from SEQSIM)')
        k = write_contacts(f'{name}_SEQSIM.contacts', sim, start)
        print(f'{name}_SEQSIM.contacts: {k} pairs')


if __name__ == '__main__':
    main()
