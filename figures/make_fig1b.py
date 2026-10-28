"""make_fig1b.py -- Fig 1B at 2 kb bins (matching SEQSIM). Run next to fig1b_microc_H1.npz.
Writes fig1b_microc_map.png (one pixel per bin) and f1b_svg.frag (overlay for the Fig 1 board)."""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
R = 2000
d = np.load('fig1b_microc_H1.npz')
B = d[f'bal_{R}'].astype(float); raw = d[f'raw_{R}'].astype(float); w = d[f'w_{R}']; s = d[f'starts_{R}']
W0, W1 = 11_800_000, 14_400_000; L0, L1 = 12_771_346, 13_423_345; U0, U1 = 12_954_411, 13_423_345
with np.errstate(divide='ignore', invalid='ignore'):
    L = np.log10(B)
fin = np.isfinite(L) & (B > 0)
lo, hi = np.percentile(L[fin], [5, 99.3])
cm = LinearSegmentedColormap.from_list('jb', ['#ffffff', '#f6b3a3', '#e0402b', '#8e1016'])
rgb = cm(np.clip(np.nan_to_num((L - lo) / (hi - lo), nan=0), 0, 1))[:, :, :3]
grey = np.array([0.914, 0.898, 0.855])
rgb[~np.isfinite(L)] = grey; rgb[B == 0] = [1, 1, 1]
nanrow = np.isnan(w); rgb[nanrow, :] = grey; rgb[:, nanrow] = grey
Image.fromarray((rgb * 255).astype(np.uint8)).save('fig1b_microc_map.png', optimize=True)
MX, MY, MS = 74, 246, 216
x = lambda g: MX + (g - W0) / (W1 - W0) * MS
o = []
o.append(f'<image href="@HREF@" x="{MX}" y="{MY}" width="{MS}" height="{MS}" preserveAspectRatio="none" style="image-rendering:pixelated"></image>')
o.append(f'<rect x="{MX}" y="{MY}" width="{MS}" height="{MS}" fill="none" stroke="#d9d4c7"></rect>')
o.append(f'<rect x="{x(L0):.1f}" y="{MY+(L0-W0)/(W1-W0)*MS:.1f}" width="{(L1-L0)/(W1-W0)*MS:.1f}" height="{(L1-L0)/(W1-W0)*MS:.1f}" fill="none" stroke="#23231f" stroke-width="1.3" stroke-dasharray="4 2"><title>PRAMEF/HNRNPCL locus, chr1:12,771,346–13,423,345</title></rect>')
for mb in (12.0, 12.5, 13.0, 13.5, 14.0):
    g = mb * 1e6; xx = x(g); yy = MY + (g - W0) / (W1 - W0) * MS
    o.append(f'<line x1="{xx:.1f}" y1="{MY+MS}" x2="{xx:.1f}" y2="{MY+MS+3}" stroke="#8a8578"></line><text x="{xx:.1f}" y="{MY+MS+12}" text-anchor="middle" font-size="8.5" fill="#6b675c">{mb:.1f}</text>')
    o.append(f'<line x1="{MX-3}" y1="{yy:.1f}" x2="{MX}" y2="{yy:.1f}" stroke="#8a8578"></line><text x="{MX-5}" y="{yy+3:.1f}" text-anchor="end" font-size="8.5" fill="#6b675c">{mb:.1f}</text>')
o.append(f'<text x="{MX+MS+4}" y="{MY+MS+12}" font-size="8.5" fill="#6b675c">Mb</text>')
# coverage track: raw read pairs per 2 kb bin, log scale; two paths (kept / discarded bins)
TY, TH = 206, 28
cov = raw.sum(1); lc = np.log10(np.maximum(cov, 1)); top = np.log10(cov.max()); bw = MS / len(s)
paths = {False: [], True: []}
for i, c in enumerate(lc):
    h = TH * c / top
    if h <= 0: continue
    paths[bool(nanrow[i])].append(f'M{MX+i*bw:.2f} {TY+TH:.2f}v{-h:.2f}h{bw+0.05:.2f}v{h:.2f}z')
o.append(f'<path d="{"".join(paths[False])}" fill="#6b675c"><title>read pairs per 2 kb bin (log scale)</title></path>')
o.append(f'<path d="{"".join(paths[True])}" fill="#b8b3a6"><title>bins discarded by balancing</title></path>')
o.append(f'<line x1="{MX}" y1="{TY+TH}" x2="{MX+MS}" y2="{TY+TH}" stroke="#8a8578" stroke-width="0.6"></line>')
o.append(f'<text x="{MX-5}" y="{TY+12}" text-anchor="end" font-size="8.5" fill="#6b675c">reads</text><text x="{MX-5}" y="{TY+22}" text-anchor="end" font-size="8.5" fill="#6b675c">(log)</text>')
o.append(f'<rect x="{x(U0):.1f}" y="{TY+TH+3}" width="{x(U1)-x(U0):.1f}" height="4" rx="1" fill="#b0561a"><title>no unique mapping (panel A): 12,954,411–13,423,345</title></rect>')
o.append(f'<rect x="{MX+MS-104}" y="{TY-11}" width="9" height="7" fill="#e9e5da" stroke="#b8b3a6" stroke-width="0.6"></rect><text x="{MX+MS-92}" y="{TY-5}" font-size="8.5" fill="#6b675c">bins discarded</text>')
open('f1b_svg.frag', 'w').write('\n'.join(o))
inl = (s >= L0 - R) & (s < L1)
print(f'{R} bp bins: locus discarded {nanrow[inl].sum()}/{inl.sum()} ({nanrow[inl].mean()*100:.0f}%), flanks discarded {nanrow[~inl].sum()}/{(~inl).sum()}')
print(f'median read pairs per bin: locus {np.median(cov[inl]):.0f}, flanks {np.median(cov[~inl]):.0f}, ratio {np.median(cov[~inl])/np.median(cov[inl]):.1f}')
print('contrast lo hi', lo, hi)
