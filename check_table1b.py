"""check_table1b.py -- re-export the ChIP evidence of Table 1 (sites C1-C9) directly from the UCSC API.

The ENCODE cell-type counts and the ReMap factors of Table 1 were first transcribed by hand.
This script fetches them again for every site and writes a table
you can compare cell by cell with Table 1.

Run in MaaD-TAD (needs requests; internet):
    python check_table1b.py
Writes table1_recheck.tsv and prints the same table.

For each site (1-based GRCh38 interval, as in Table 1) it reports
  ENCODE (encRegTfbsClustered): RAD21, SMC3 and CTCF clusters overlapping the interval, with the
      number of cell types (sourceCount); if several clusters overlap, the largest count is shown
  ReMap 2022: cohesin subunits (RAD21, SMC3, SMC1A, STAG1, STAG2) and CTCF with a peak overlapping
  ENCODE cCREs (encodeCcreCombined) overlapping the interval, with their class
Columns *_300 repeat the ENCODE counts within 300 bp of the interval (the distance used for 'predicted' sites).
"""
import csv, re, sys, time
import requests

API = 'https://api.genome.ucsc.edu'
SITES = [('C1', 12779048, 12779264), ('C2', 12781561, 12781777), ('C3', 12783736, 12784186),
         ('C4', 12806339, 12806789), ('C5', 13228393, 13228547), ('C6', 13277160, 13277305),
         ('C7', 13356712, 13357255), ('C8', 13368533, 13368929), ('C9', 13422468, 13422998)]
COHESIN = ('RAD21', 'SMC3', 'SMC1A', 'STAG1', 'STAG2')
LO, HI = 12770000, 13425000


def get(path, **params):
    for i in range(3):
        try:
            r = requests.get(API + path, params=params, timeout=300)
            r.raise_for_status(); return r.json()
        except Exception as e:
            if i == 2: print('  failed:', path, params.get('track'), e); return None
            time.sleep(3)


def track(name):
    j = get('/getData/track', genome='hg38', track=name, chrom='chr1', start=LO, end=HI, maxItemsOutput=1000000)
    if not j: return []
    rows = j.get(name, [])
    return rows.get('chr1', []) if isinstance(rows, dict) else rows


def remap_track():
    j = get('/list/tracks', genome='hg38', trackLeavesOnly=1)
    names = [k for k in (j or {}).get('hg38', {}) if 'remap' in k.lower()]
    print('ReMap tracks on hg38:', names or 'none found')
    for pref in ('reMap2022', 'ReMap2022', 'remap2022'):
        if pref in names: return pref
    return names[0] if names else None


def overlap(r, s, e, pad=0):
    a = int(r.get('chromStart', r.get('start', 0))) + 1; b = int(r.get('chromEnd', r.get('end', 0)))
    return b >= s - pad and a <= e + pad


def main():
    tf = track('encRegTfbsClustered'); print('encRegTfbsClustered rows:', len(tf))
    cc = track('encodeCcreCombined'); print('encodeCcreCombined rows:', len(cc))
    rn = remap_track(); rm = track(rn) if rn else []; print('ReMap rows:', len(rm))
    for r in rm[:3]: print('  ReMap row example:', {k: r[k] for k in list(r)[:8]})
    out = []
    for sid, s, e in SITES:
        row = {'site': sid, 'interval': f'{s:,}-{e:,}'}
        for pad, tag in ((0, ''), (300, '_300')):
            for f in ('RAD21', 'SMC3', 'CTCF'):
                hits = [int(r.get('sourceCount', 0)) for r in tf if r.get('name') == f and overlap(r, s, e, pad)]
                row[f'ENCODE_{f}{tag}'] = max(hits) if hits else ''
        fac = set()
        for r in rm:
            if overlap(r, s, e):
                txt = ' '.join(str(v) for k, v in r.items() if k not in ('chrom', 'chromStart', 'chromEnd'))
                toks = set(re.split(r'[^A-Za-z0-9]+', txt.upper()))
                for f in COHESIN + ('CTCF',):
                    if f in toks: fac.add(f)
        row['ReMap_cohesin'] = ', '.join(f for f in COHESIN if f in fac) or '-'
        row['ReMap_CTCF'] = 'yes' if 'CTCF' in fac else '-'
        row['cCRE_300'] = '; '.join(str(r.get('name')) for r in cc if overlap(r, s, e, 300)) or '-'
        row['cCRE'] = '; '.join(f"{r.get('name')} ({r.get('ucscLabel') or r.get('encodeLabel') or r.get('ccre', '')})"
                                for r in cc if overlap(r, s, e)) or '-'
        out.append(row)
    cols = list(out[0])
    with open('table1_recheck.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t'); w.writeheader(); w.writerows(out)
    print('\n' + '\t'.join(cols))
    for r in out: print('\t'.join(str(r[c]) for c in cols))
    print('\nwrote table1_recheck.tsv  (compare with Table 1: ENCODE counts in parentheses, ReMap factors, cCRE accessions)')


if __name__ == '__main__':
    main()
