import json, csv, math, re
S,E_=12771346,13423345
cx,cy,R=548,332,102
T0,SPAN=18.0,324.0
ang=lambda g: math.radians(T0+(g-S)/(E_-S)*SPAN)
P=lambda g,r:(cx+r*math.sin(ang(g)), cy-r*math.cos(ang(g)))
f=lambda p:f'{p[0]:.1f} {p[1]:.1f}'
d=json.load(open('data/homology_restraints.json'))
G=[r for r in csv.reader(open('data/genes.tsv'),delimiter='\t')]
GREEN,MAG,INV,DIR,INK,MUT='#2f8f5b','#a8418a','#d9731f','#3f6fa8','#23231f','#6b675c'
o=[]
def arcpath(a,b,r,sweep=1):
    large=1 if abs(math.degrees(ang(b)-ang(a)))>180 else 0
    return f'A {r} {r} 0 {large} {sweep} {f(P(b,r))}'
# ribbons (draw low identity first)
rr=R-4
for b in sorted(d['blocks'],key=lambda b:b['ident']):
    a0,a1=b['armA']; c0,c1=b['armC']
    if a0>c0: a0,a1,c0,c1=c0,c1,a0,a1
    col=INV if b['type']=='inverted' else DIR
    op=0.30+0.45*(b['ident']-97.9)/(100-97.9)
    C=f'{cx} {cy}'
    if b['type']=='inverted':   # untwisted: a0->a1, a1->c0, c0->c1, c1->a0 (reverse correspondence)
        path=f'M {f(P(a0,rr))} {arcpath(a0,a1,rr)} Q {C} {f(P(c0,rr))} {arcpath(c0,c1,rr)} Q {C} {f(P(a0,rr))} Z'
    else:                       # twisted: a1->c1, c0->a0 (same orientation)
        path=f'M {f(P(a0,rr))} {arcpath(a0,a1,rr)} Q {C} {f(P(c1,rr))} {arcpath(c1,c0,rr,0)} Q {C} {f(P(a0,rr))} Z'
    o.append(f'<path d="{path}" fill="{col}" fill-opacity="{op:.2f}" stroke="{col}" stroke-width="0.6"><title>{b["id"]} · {b["type"]} · {b["ident"]}% identity · arms {b["armA"][0]:,}–{b["armA"][1]:,} and {b["armC"][0]:,}–{b["armC"][1]:,}</title></path>')
# backbone ring
o.append(f'<path d="M {f(P(S,R))} {arcpath(S,E_,R)}" fill="none" stroke="#8a8578" stroke-width="1.2"></path>')
# genes as ring arcs
for name,bt,a,b_,st,tx in G:
    a,b_=int(a),int(b_)
    if bt=='protein_coding':
        col=MAG if name.startswith('HNRNPCL') else GREEN
        o.append(f'<path d="M {f(P(a,R+5))} {arcpath(a,b_,R+5)}" fill="none" stroke="{col}" stroke-width="9"><title>{name} ({"+" if int(st)>0 else "−"})</title></path>')
    else:
        o.append(f'<path d="M {f(P(a,R+5))} {arcpath(a,b_,R+5)}" fill="none" stroke="#b8b3a6" stroke-width="5"><title>{name} (pseudogene)</title></path>')
# ticks every 100 kb + ends
for mb in [12.8,12.9,13.0,13.1,13.2,13.3,13.4]:
    g=mb*1e6; p1,p2,pl=P(g,R+11),P(g,R+15),P(g,R+36)
    o.append(f'<line x1="{p1[0]:.1f}" y1="{p1[1]:.1f}" x2="{p2[0]:.1f}" y2="{p2[1]:.1f}" stroke="#8a8578"></line><text x="{pl[0]:.1f}" y="{pl[1]+3:.1f}" text-anchor="middle" font-size="9.5" fill="{MUT}">{mb:.1f}</text>')
for lab,deg,anc in (('start →',6,'start'),('← end',354,'end')):
    th=math.radians(deg); x_=cx+(R+30)*math.sin(th); y_=cy-(R+30)*math.cos(th)
    o.append(f'<text x="{x_:.1f}" y="{y_+3:.1f}" text-anchor="{anc}" font-size="9.5" font-style="italic" fill="{INK}">{lab}</text>')
# assembly gap (all-N) in GRCh38
GA,GB=12954385,13004384
o.append(f'<path d="M {f(P(GA,R))} {arcpath(GA,GB,R)}" fill="none" stroke="#ffffff" stroke-width="4"></path>')
o.append(f'<path d="M {f(P(GA,R))} {arcpath(GA,GB,R)}" fill="none" stroke="#8a8578" stroke-width="1.2" stroke-dasharray="2 2"><title>assembly gap (50 kb of N) {GA:,}–{GB:,}</title></path>')
pg=P((GA+GB)//2,R-12); o.append(f'<text x="{pg[0]:.1f}" y="{pg[1]+3:.1f}" text-anchor="middle" font-size="9" font-style="italic" fill="{MUT}">gap</text>')
# sites: ChIP-supported cohesin (filled triangle), CTCF only (open triangle), predicted mirror copies (diamonds, new colour)
NEWC='#c1272d'
coh=[(12779156,'RAD21'),(12781669,'RAD21'),(12783961,'RAD21 (T_I2)'),(12806564,'RAD21'),(13228470,'SMC3, weak'),
     (13277233,'RAD21 (T_I5)'),(13356984,'RAD21+SMC3 (F5)'),(13368731,'RAD21 (T_I7)'),(13422733,'RAD21+SMC3')]
ctc=[12819268,12823199,12863287,12921867,12924608,13005906,13247275,13279908,13324390,13355034]
new=json.load(open('data/new_sites.json'))
def tri(g,col,fill,title):
    tip=P(g,R+10.5); bl=P(g-2600,R+19); br=P(g+2600,R+19)
    return f'<path d="M {f(tip)} L {f(bl)} L {f(br)} Z" fill="{fill}" stroke="{col}" stroke-width="1"><title>{title}</title></path>'
def dia(g,fill,title):
    c_=P(g,R+24); r_=4.6
    pts=[(c_[0],c_[1]-r_),(c_[0]+r_,c_[1]),(c_[0],c_[1]+r_),(c_[0]-r_,c_[1])]
    return f'<path d="M {" L ".join(f(q) for q in pts)} Z" fill="{fill}" stroke="{NEWC}" stroke-width="1.2"><title>{title}</title></path>'
for g in ctc: o.append(tri(g,INK,'#ffffff',f'CTCF ChIP peak, no cohesin ~{g:,}'))
for g,t in coh: o.append(tri(g,INK,INK,f'cohesin ChIP peak: {t} ~{g:,}'))
for v in new:
    fill=NEWC if v['kind']=='cohesin family' else '#ffffff'
    o.append(dia(v['pos'],fill,f"predicted mirror site {v['pos']:,}: copy of {','.join(v['sources'])} ({v['kind']}), {v['ident']*100:.1f}% identity, CTCF motif score {v['motif']['rel']}, no ChIP peak"))
# legend (right)
lx,ly=716,226
def rib(x,y,col,tw):
    if tw: return f'<path d="M {x} {y-6} L {x+22} {y+4} L {x+22} {y-6} L {x} {y+4} Z" fill="{col}" fill-opacity="0.6" stroke="{col}" stroke-width="0.6"></path>'
    return f'<rect x="{x}" y="{y-6}" width="22" height="10" fill="{col}" fill-opacity="0.6" stroke="{col}" stroke-width="0.6"></rect>'
o.append(rib(lx,ly,INV,False)+f'<text x="{lx+30}" y="{ly+2}" font-size="10.5" fill="{INK}">inverted repeat (untwisted)</text>')
o.append(rib(lx,ly+20,DIR,True)+f'<text x="{lx+30}" y="{ly+22}" font-size="10.5" fill="{INK}">direct repeat (twisted)</text>')
o.append(f'<text x="{lx}" y="{ly+40}" font-size="9.5" fill="{MUT}">darker ribbon = higher identity (97.9 → 99.98%)</text>')
y=ly+62
o.append(f'<rect x="{lx}" y="{y-7}" width="10" height="8" fill="{GREEN}"></rect><text x="{lx+15}" y="{y}" font-size="10" fill="{INK}">PRAMEF</text><rect x="{lx+66}" y="{y-7}" width="10" height="8" fill="{MAG}"></rect><text x="{lx+81}" y="{y}" font-size="10" fill="{INK}">HNRNPCL</text><rect x="{lx+136}" y="{y-6}" width="10" height="6" fill="#b8b3a6"></rect><text x="{lx+151}" y="{y}" font-size="10" fill="{INK}">pseudo</text>')
y+=22
o.append(f'<text x="{lx}" y="{y}" font-size="10" font-weight="600" fill="{INK}">ChIP-seq peaks (ENCODE, ReMap)</text>')
y+=16; o.append(f'<path d="M {lx+5} {y} L {lx} {y-8} L {lx+10} {y-8} Z" fill="{INK}"></path><text x="{lx+15}" y="{y}" font-size="10" fill="{INK}">cohesin + CTCF (9)</text><path d="M {lx+120} {y} L {lx+115} {y-8} L {lx+125} {y-8} Z" fill="#fff" stroke="{INK}"></path><text x="{lx+130}" y="{y}" font-size="10" fill="{INK}">CTCF only (10)</text>')
y+=22
o.append(f'<text x="{lx}" y="{y}" font-size="10" font-weight="600" fill="{NEWC}">Predicted mirror sites, no peak (15)</text>')
y+=16; o.append(f'<path d="M {lx+5} {y-9} L {lx+10} {y-4} L {lx+5} {y+1} L {lx} {y-4} Z" fill="{NEWC}" stroke="{NEWC}"></path><text x="{lx+15}" y="{y}" font-size="10" fill="{INK}">copy of a cohesin site (6)</text>')
y+=15; o.append(f'<path d="M {lx+5} {y-9} L {lx+10} {y-4} L {lx+5} {y+1} L {lx} {y-4} Z" fill="#fff" stroke="{NEWC}" stroke-width="1.2"></path><text x="{lx+15}" y="{y}" font-size="10" fill="{INK}">copy of a CTCF-only site (9)</text>')
y+=15; o.append(f'<text x="{lx}" y="{y}" font-size="9.5" fill="{MUT}">≥94% identical, CTCF motif intact;</text>')
y+=12; o.append(f'<text x="{lx}" y="{y}" font-size="9.5" fill="{MUT}">14 of 15 lie in unmappable sequence</text>')
newC='\n'.join(o)
open('fig1c.svg','w').write('<svg xmlns="http://www.w3.org/2000/svg" width="928" height="600" viewBox="0 0 928 600"><g font-family="IBM Plex Sans, sans-serif" font-size="12" fill="#23231f">'+newC+'</g></svg>')
print('wrote fig1c.svg')
