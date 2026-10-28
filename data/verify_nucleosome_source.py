"""verify_nucleosome_source.py -- which UCSC nucleosome track did our 490 measured dyads come from?

Run in MaaD-TAD\\data (next to nucleosomes_GRCh38_resolved.bed):
    pip install pyliftover requests
    python verify_nucleosome_source.py

What it does
  1. Reads the 490 'real' dyads (GRCh38) from nucleosomes_GRCh38_resolved.bed.
  2. Lists UCSC nucleosome tracks for hg19 and hg38 (UCSC REST API) whose name or label mentions
     nucleosomes / MNase, keeping GM12878 ones first.
  3. For hg19 tracks: lifts the dyads GRCh38 -> hg19 (pyliftover, UCSC chain) and reads the track
     over the lifted window. For hg38 tracks: reads the track directly (no liftover).
  4. For each track reports how well our dyads match it:
       - share of dyads within 10 bp of a local maximum of the signal (random positions as control)
       - Spearman correlation between our 'score' column and the track value at the dyad
       - share of dyads whose score equals the track value at the dyad (an exact match = this is the source)
Writes nucleosome_source_check.txt next to the bed and prints the same report.
"""
import json, random, sys
import requests
from pyliftover import LiftOver

BED = 'nucleosomes_GRCh38_resolved.bed'
API = 'https://api.genome.ucsc.edu'
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

# 1. our measured dyads
dy = []
for line in open(BED):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if f[7] == 'real':
        dy.append((int(f[6]), float(f[4])))
say(f'{len(dy)} measured dyads in {BED} (GRCh38)')

# 2. candidate tracks
def tracks(genome):
    r = requests.get(f'{API}/list/tracks', params={'genome': genome, 'trackLeavesOnly': 1}, timeout=120); r.raise_for_status()
    d = r.json()[genome]; hits = []
    for k, v in d.items():
        lab = (k + ' ' + str(v.get('shortLabel', '')) + ' ' + str(v.get('longLabel', ''))).lower()
        if ('nsome' in lab or 'nucleosom' in lab or 'mnase' in lab) and str(v.get('type', '')).split()[0] in ('bigWig', 'wig', 'bedGraph', 'bigBed', 'bed', 'narrowPeak', 'broadPeak'):
            hits.append((0 if 'gm12878' in lab else 1, k, v.get('longLabel', ''), v.get('type', '')))
    return sorted(hits)

def get(genome, track, chrom, start, end):
    r = requests.get(f'{API}/getData/track', params={'genome': genome, 'track': track, 'chrom': chrom,
                     'start': start, 'end': end, 'maxItemsOutput': 1000000}, timeout=300)
    if r.status_code != 200: return None
    j = r.json(); rows = j.get(track) or j.get(chrom) or []
    if isinstance(rows, dict): rows = rows.get(chrom, [])
    return rows

def per_bp(rows, start, end):
    """value at each base of [start, end) from bigWig/bedGraph-style rows (start, end, value)."""
    import numpy as np
    v = np.full(end - start, np.nan)
    for r in rows:
        s, e = int(r.get('start', r.get('chromStart', 0))), int(r.get('end', r.get('chromEnd', 0)))
        val = r.get('value', r.get('score', r.get('signalValue')))
        if val is None: continue
        a, b = max(s, start) - start, min(e, end) - start
        if b > a: v[a:b] = float(val)
    return v

def score(genome, track, pos, scores):
    import numpy as np
    from scipy.stats import spearmanr
    lo, hi = min(pos) - 5000, max(pos) + 5000
    rows = get(genome, track, 'chr1', lo, hi)
    if not rows: return None
    v = per_bp(rows, lo, hi)
    if np.isnan(v).all(): return None
    w = np.nan_to_num(v, nan=-1)
    def near_peak(p):
        i = p - lo
        if i < 75 or i >= len(w) - 75: return False
        win = w[i - 75:i + 76]; j = int(np.argmax(win)) - 75
        return abs(j) <= 10 and win.max() > 0
    hit = np.mean([near_peak(p) for p in pos])
    rnd = random.Random(1); ctrl = np.mean([near_peak(rnd.randrange(lo + 100, hi - 100)) for _ in range(2000)])
    at = np.array([v[p - lo] for p in pos])
    ok = np.isfinite(at)
    rho = spearmanr(np.array(scores)[ok], at[ok]).correlation if ok.sum() > 10 else float('nan')
    exact = np.mean(np.isclose(np.array(scores)[ok], at[ok], atol=1e-3)) if ok.any() else 0
    return hit, ctrl, rho, exact, ok.sum()

pos38 = [p for p, _ in dy]; sc = [s for _, s in dy]
# 3a. hg38 tracks (no liftover)
say('\n== hg38 tracks (positions used as is) ==')
for _, k, lab, typ in tracks('hg38')[:12]:
    r = score('hg38', k, pos38, sc)
    say(f'{k} [{typ}] {lab}\n    ' + ('no data over the locus' if r is None else
        f'dyads at a local signal max (+-10 bp): {r[0]:.0%} (random {r[1]:.0%}) | Spearman score vs signal {r[2]:.2f} | exact value match {r[3]:.0%} (n={r[4]})'))

# 3b. hg19 tracks (lift GRCh38 -> hg19 first)
lo = LiftOver('hg38', 'hg19')
lifted = [(lo.convert_coordinate('chr1', p - 1), s) for p, s in dy]
keep = [(c[0][1] + 1, s) for c, s in lifted if c and c[0][0] == 'chr1']
say(f'\n== hg19 tracks: {len(keep)} of {len(dy)} dyads lift to hg19 chr1 '
    f'({min(p for p, _ in keep):,}-{max(p for p, _ in keep):,}) ==')
pos19 = [p for p, _ in keep]; sc19 = [s for _, s in keep]
for _, k, lab, typ in tracks('hg19')[:12]:
    r = score('hg19', k, pos19, sc19)
    say(f'{k} [{typ}] {lab}\n    ' + ('no data over the locus' if r is None else
        f'dyads at a local signal max (+-10 bp): {r[0]:.0%} (random {r[1]:.0%}) | Spearman score vs signal {r[2]:.2f} | exact value match {r[3]:.0%} (n={r[4]})'))

say('\nHow to read it: the source track is the one with most dyads at a local maximum (far above random),'
    ' a high Spearman, and ideally exact value matches. If it is an hg19 track, the liftOver step is confirmed.')
open('nucleosome_source_check.txt', 'w').write('\n'.join(out) + '\n')
print('\nwrote nucleosome_source_check.txt')
