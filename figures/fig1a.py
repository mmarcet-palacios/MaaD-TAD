import csv, re, json, numpy as np
S,E_=12771346,13423345
X0,X1=104,912; W=X1-X0
x=lambda g: X0+(g-S)/(E_-S)*W
G=[r for r in csv.reader(open('data/genes.tsv'),delimiter='\t')]
EX={}
for t,a,b in csv.reader(open('data/exons.tsv'),delimiter='\t'): EX.setdefault(t,[]).append((int(a),int(b)))
GREEN,MAG,GREY,BLUE,RUST='#2f8f5b','#a8418a','#b8b3a6','#3f6fa8','#b0561a'
INK,MUTED='#23231f','#6b675c'
short=lambda nm: nm.replace('PRAMEF','').replace('HNRNPCL','CL')
out=[]; T='font-size="18" font-weight="600"'; M=f'fill="{MUTED}"'
out.append(f'<text x="12" y="24" {T}>A</text><text x="40" y="24" {M}>genes (Ensembl, GRCh38) over the nucleosome data, to scale</text>')
out.append(f'<line x1="{X0}" y1="48" x2="{X1}" y2="48" stroke="#8a8578"></line><text x="96" y="40" text-anchor="end" font-size="10" {M}>Mb</text>')
for mb in np.arange(12.8,13.45,0.1):
    xx=x(mb*1e6); out.append(f'<line x1="{xx:.1f}" y1="44" x2="{xx:.1f}" y2="48" stroke="#8a8578"></line><text x="{xx:.1f}" y="40" text-anchor="middle" font-size="10" {M}>{mb:.1f}</text>')
rows={1:80,-1:114}
out.append(f'<text x="{X0-8}" y="{rows[1]+4}" text-anchor="end" font-size="10.5" {M}>+ strand</text><text x="{X0-8}" y="{rows[-1]+4}" text-anchor="end" font-size="10.5" {M}>− strand</text>')
for s,y in rows.items(): out.append(f'<line x1="{X0}" y1="{y}" x2="{X1}" y2="{y}" stroke="#eee9dd"></line>')
labels={1:[],-1:[]}
for name,bt,a,b,st,tx in sorted(G,key=lambda r:int(r[2])):
    a,b,st=int(a),int(b),int(st); y=rows[st]
    if bt=='protein_coding':
        col=MAG if name.startswith('HNRNPCL') else GREEN
        g=[f'<g><title>{name} · {a:,}–{b:,} ({"+" if st>0 else "−"}) · {tx}</title>',
           f'<line x1="{x(a):.1f}" y1="{y}" x2="{x(b):.1f}" y2="{y}" stroke="{col}" stroke-width="1.5"></line>']
        for ea,eb in EX[tx]:
            g.append(f'<rect x="{x(ea):.2f}" y="{y-6}" width="{max(1.2,x(eb)-x(ea)):.2f}" height="12" rx="1" fill="{col}"></rect>')
        tip=x(b)+4 if st>0 else x(a)-4
        g.append(f'<path d="M {tip-3*st:.1f} {y-4} L {tip:.1f} {y} L {tip-3*st:.1f} {y+4}" fill="none" stroke="{col}" stroke-width="1.3"></path></g>')
        out+=g; labels[st].append(((x(a)+x(b))/2,short(name),INK))
    else:
        out.append(f'<rect x="{x(a):.2f}" y="{y-4}" width="{x(b)-x(a):.2f}" height="8" rx="1" fill="#ffffff" stroke="#8a8578" stroke-width="1.2"><title>{name} ({bt.replace("_"," ")}) · {a:,}–{b:,}</title></rect>')
        labels[st].append(((x(a)+x(b))/2,short(name),'#8a8578'))
for st,L in labels.items():
    y0=rows[st]-11 if st>0 else rows[st]+20
    last=[-1e9,-1e9,-1e9]
    for cx,t,col in sorted(L):
        w=6.0*len(t)
        tier=next((k for k in range(3) if cx-w/2>last[k]+2),2)
        last[tier]=cx+w/2
        yy=y0-(9*tier if st>0 else -9*tier)
        it=' font-style="italic"' if t.endswith('P') else ''
        out.append(f'<text x="{cx:.1f}" y="{yy}" text-anchor="middle" font-size="10" fill="{col}"{it}>{t}</text>')
# nucleosome tracks
R=[l.rstrip('\n').split('\t') for l in open('data/nucleosomes_GRCh38_resolved.bed') if not l.startswith('#')]
R.sort(key=lambda r:int(r[6]))
n=160; real=[0]*n; tot=[0]*n
for r in R:
    c=int(r[6]); b=min(n-1,max(0,int((c-S)/(E_-S)*n))); tot[b]+=1; real[b]+=r[7]=='real'
bw=W/n; YN=158
out.append(f'<text x="{X0-8}" y="{YN+15}" text-anchor="end" font-size="10.5" {M}>measured</text><text x="{X0-8}" y="{YN+27}" text-anchor="end" font-size="10.5" {M}>dyads</text>')
out.append(f'<rect x="{X0}" y="{YN}" width="{W}" height="30" fill="#f4f2ec"></rect>')
for i in range(n):
    f=real[i]/max(1,tot[i])
    if f>0: out.append(f'<rect x="{X0+i*bw:.2f}" y="{YN+30-30*f:.2f}" width="{bw-0.6:.2f}" height="{30*f:.2f}" fill="{BLUE}"><title>{(S+i*(E_-S)/n)/1e6:.3f} Mb: {real[i]} of {tot[i]} nucleosomes measured</title></rect>')
# unmappable runs (exact, from the nucleosome records)
runs=[]; cur=None
for r in R:
    u=r[8]=='unmapped_assembly'; s0,e0=int(r[1]),int(r[2])
    if u and cur: cur[1]=e0; cur[2]+=1
    elif u: cur=[s0,e0,1]
    elif cur: runs.append(cur); cur=None
if cur: runs.append(cur)
nu=sum(c for _,_,c in runs); ntot=len(R); bp=sum(e-s for s,e,_ in runs)
out.append(f'<text x="{X0-8}" y="{YN+48}" text-anchor="end" font-size="10.5" {M}>unmappable</text>')
for s0,e0,c in runs:
    out.append(f'<rect x="{x(s0):.2f}" y="{YN+38}" width="{x(e0)-x(s0):.2f}" height="14" rx="2" fill="{RUST}"><title>{s0:,}–{e0:,}: {c:,} nucleosomes</title></rect>')
big=max(runs,key=lambda r:r[1]-r[0])
out.append(f'<text x="{x(big[0])+8:.1f}" y="{YN+48.5}" font-size="10.5" font-weight="600" fill="#ffffff">{nu/ntot:.0%} of nucleosomes ({nu:,} of {ntot:,}) · {bp/1000:.0f} kb, {bp/(E_-S+1):.0%} of the locus</text>')
# legend: colour boxes and their meaning only (details live in the figure caption)
y=YN+76
def swatch(xx,yy,col,kind='box'):
    return (f'<rect x="{xx}" y="{yy-8}" width="12" height="9" rx="1" fill="{col}"></rect>' if kind=='box'
            else f'<rect x="{xx}" y="{yy-7}" width="12" height="7" rx="1" fill="#fff" stroke="#8a8578" stroke-width="1.2"></rect>')
items=[(GREEN,'box','PRAMEF gene'),(MAG,'box','HNRNPCL gene'),('#8a8578','outline','pseudogene'),
       (BLUE,'box','measured dyads (fraction per 4 kb)'),(RUST,'box','no unique mapping')]
xx=40
for col,kind,txt in items:
    out.append(swatch(xx,y,col,kind)+f'<text x="{xx+17}" y="{y}" font-size="10.5" {M}>{txt}</text>')
    xx+=17+6.0*len(txt)+22
end=y+14
print('panel A ends at',end,'| unmappable',nu,ntot,bp)
panelA='\n'.join(out)

open('fig1a.svg','w').write('<svg xmlns="http://www.w3.org/2000/svg" width="928" height="600" viewBox="0 0 928 600"><g font-family="IBM Plex Sans, sans-serif" font-size="12" fill="#23231f">'+panelA+'</g></svg>')
print('wrote fig1a.svg')
