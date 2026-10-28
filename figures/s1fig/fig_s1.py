import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
a = np.load(os.path.join(HERE, 'anch_s1fig.npz')); f = np.load(os.path.join(HERE, 'free_s1fig.npz'))
v = [-0.529923, 0.789860, -0.308712, -0.606986, -0.607489, -0.512372, -0.592241, -0.084134, 0.801356]
A = np.array(v).reshape(3, 3); o = a['origins']; c = (o - o.mean(0)) @ A; X, Y = c[:, 0], -c[:, 1]
real = a['src'] == 'real'
out = []; E = out.append
# ---------- Panel A ----------
E('<text x="12" y="28" font-size="18" font-weight="600">A</text><text x="40" y="28" fill="#6b675c">one point per nucleosome (3,359)</text>')
x0, y0, W, H = 24, 44, 360, 400
s = min(W / (X.max() - X.min()), H / (Y.max() - Y.min()))
px = x0 + (W - (X.max() - X.min()) * s) / 2 + (X - X.min()) * s
py = y0 + (H - (Y.max() - Y.min()) * s) / 2 + (Y - Y.min()) * s
E('<g fill="#a39e90" fill-opacity="0.75">' + ''.join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.5"/>' for x, y in zip(px[~real], py[~real])) + '</g>')
E('<g fill="#3f6fa8" stroke="#ffffff" stroke-width="0.5">' + ''.join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.1"/>' for x, y in zip(px[real], py[real])) + '</g>')
bar = 1000 * s  # 100 nm = 1000 A
E(f'<line x1="30" y1="462" x2="{30+bar:.1f}" y2="462" stroke="#23231f" stroke-width="2"/><text x="{36+bar:.1f}" y="466" font-size="9.5" fill="#6b675c">100 nm</text>')
E(f'<circle cx="34" cy="484" r="3" fill="#3f6fa8"/><text x="42" y="487.5" font-size="10">measured dyad ({real.sum()})</text>')
E(f'<circle cx="170" cy="484" r="3" fill="#a39e90"/><text x="178" y="487.5" font-size="10">placeholder ({(~real).sum():,})</text>')
E('<line x1="404" y1="40" x2="404" y2="500" stroke="#e4dfd2"/>')
# ---------- Panel B ----------
corr = a['corr']; med = np.median(corr); p95 = np.percentile(corr, 95)
E('<text x="416" y="28" font-size="18" font-weight="600">B</text><text x="444" y="28" fill="#6b675c">steering rotation per nucleosome</text>')
bx0, bx1, by0, by1 = 470, 900, 50, 200
edges = np.arange(0, 20.5, 0.5); h, _ = np.histogram(corr[corr < 20], edges); over = int((corr >= 20).sum())
ymax = 800
XB = lambda d: bx0 + (d / 20) * (bx1 - bx0)
YB = lambda n: by1 - (n / ymax) * (by1 - by0)
for t in (0, 200, 400, 600, 800):
    E(f'<line x1="{bx0-3}" y1="{YB(t):.1f}" x2="{bx0}" y2="{YB(t):.1f}" stroke="#8a8578"/><text x="{bx0-6}" y="{YB(t)+3:.1f}" text-anchor="end" font-size="9" fill="#6b675c">{t}</text>')
    if t:
        E(f'<line x1="{bx0}" y1="{YB(t):.1f}" x2="{bx1}" y2="{YB(t):.1f}" stroke="#f0ece2"/>')
w = (bx1 - bx0) / 40
for i, n in enumerate(h):
    if n == 0:
        continue
    lo = edges[i]
    E(f'<rect x="{XB(lo)+1:.1f}" y="{YB(n):.1f}" width="{w-2:.1f}" height="{by1-YB(n):.1f}" fill="#3f6fa8"><title>{lo:.1f}–{lo+0.5:.1f}°: {n} nucleosomes</title></rect>')
E(f'<line x1="{bx0}" y1="{by1}" x2="{bx1}" y2="{by1}" stroke="#8a8578"/><line x1="{bx0}" y1="{by0}" x2="{bx0}" y2="{by1}" stroke="#8a8578"/>')
for t in (0, 5, 10, 15, 20):
    E(f'<line x1="{XB(t):.1f}" y1="{by1}" x2="{XB(t):.1f}" y2="{by1+3}" stroke="#8a8578"/><text x="{XB(t):.1f}" y="{by1+13}" text-anchor="middle" font-size="9" fill="#6b675c">{t}°</text>')
E(f'<text x="{(bx0+bx1)/2}" y="{by1+27}" text-anchor="middle" font-size="9.5" fill="#6b675c">rotation applied to follow the path (degrees per nucleosome)</text>')
E(f'<text x="436" y="{(by0+by1)/2}" font-size="9.5" fill="#6b675c" transform="rotate(-90 436 {(by0+by1)/2})" text-anchor="middle">nucleosomes</text>')
for val, lab, dy in ((med, f'median {med:.1f}°', 0), (p95, f'95th percentile {p95:.1f}°', 0)):
    E(f'<line x1="{XB(val):.1f}" y1="{by0-2}" x2="{XB(val):.1f}" y2="{by1}" stroke="#23231f" stroke-width="1" stroke-dasharray="3 3"/>')
    E(f'<text x="{XB(val)+5:.1f}" y="{by0+8+dy}" font-size="10" fill="#23231f">{lab}</text>')
E(f'<text x="{bx1}" y="{by0+40}" text-anchor="end" font-size="9.5" fill="#6b675c">{over} nucleosomes above 20° (max {corr.max():.0f}°) not shown</text>')
E('<line x1="416" y1="246" x2="916" y2="246" stroke="#e4dfd2"/>')
# ---------- Panel C ----------
E('<text x="416" y="272" font-size="18" font-weight="600">C</text><text x="444" y="272" fill="#6b675c">offset from the genomic position, along the path</text>')
cx0, cx1, cy0, cy1 = 470, 900, 312, 450
kb = (a['pos'] - 12771273) / 1000; L = 652.0; lim = 60
XC = lambda k: cx0 + k / L * (cx1 - cx0)
YC = lambda nm: cy0 + (lim - np.clip(nm, -lim, lim)) / (2 * lim) * (cy1 - cy0)
for t in (-60, -30, 0, 30, 60):
    lab = '0' if t == 0 else f'{t:+d}'.replace('-', '−')
    E(f'<line x1="{cx0-3}" y1="{YC(t):.1f}" x2="{cx0}" y2="{YC(t):.1f}" stroke="#8a8578"/><text x="{cx0-6}" y="{YC(t)+3:.1f}" text-anchor="end" font-size="9" fill="#6b675c">{lab}</text>')
    col = '#c9c4b5' if t == 0 else '#f0ece2'
    dash = ' stroke-dasharray="3 3"' if t == 0 else ''
    E(f'<line x1="{cx0}" y1="{YC(t):.1f}" x2="{cx1}" y2="{YC(t):.1f}" stroke="{col}"{dash}/>')
def poly(lag, col, extra):
    nm = lag / 10
    pts = ' '.join(f'{XC(k):.1f},{YC(v):.1f}' for k, v in zip(kb[::2], nm[::2]))
    E(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="0.9" stroke-linejoin="round"{extra}/>')
poly(f['lag_local'], '#b0561a', ' stroke-dasharray="3 1.5" stroke-opacity="0.9"')
poly(a['lag_local'], '#2f8f5b', '')
E(f'<line x1="{cx0}" y1="{cy1}" x2="{cx1}" y2="{cy1}" stroke="#8a8578"/><line x1="{cx0}" y1="{cy0}" x2="{cx0}" y2="{cy1}" stroke="#8a8578"/>')
for t in (0, 100, 200, 300, 400, 500, 600):
    E(f'<line x1="{XC(t):.1f}" y1="{cy1}" x2="{XC(t):.1f}" y2="{cy1+3}" stroke="#8a8578"/><text x="{XC(t):.1f}" y="{cy1+13}" text-anchor="middle" font-size="9" fill="#6b675c">{t}</text>')
E(f'<text x="{(cx0+cx1)/2}" y="{cy1+27}" text-anchor="middle" font-size="9.5" fill="#6b675c">position along the locus (kb)</text>')
E(f'<text x="436" y="{(cy0+cy1)/2}" font-size="9.5" fill="#6b675c" transform="rotate(-90 436 {(cy0+cy1)/2})" text-anchor="middle">offset (nm)</text>')
E(f'<path d="M {XC(649):.1f} {cy1-3} l -4 -7 h 8 z" fill="#23231f"/><text x="{XC(642):.1f}" y="{cy1-6}" text-anchor="end" font-size="9" fill="#6b675c">last nucleosomes, −140 nm (both runs; off scale)</text>')
la = np.abs(a['lag_local']) / 10; lf = np.abs(f['lag_local']) / 10
E(f'<line x1="480" y1="284" x2="498" y2="284" stroke="#2f8f5b" stroke-width="2"/><text x="504" y="287.5" font-size="10">anchored: median |offset| {np.median(la):.0f} nm, 95% within {np.percentile(la,95):.0f} nm</text>')
E(f'<line x1="480" y1="298" x2="498" y2="298" stroke="#b0561a" stroke-width="2" stroke-dasharray="3 1.5"/><text x="504" y="301.5" font-size="10">not anchored: median |offset| {np.median(lf):.0f} nm, 95% within {np.percentile(lf,95):.0f} nm</text>')
open(os.path.join(HERE, 'fig_s1_panels.svgfrag'), 'w').write('\n'.join(out))
print(len('\n'.join(out)), med, p95, over, np.median(la), np.percentile(la, 95), np.median(lf), np.percentile(lf, 95))
