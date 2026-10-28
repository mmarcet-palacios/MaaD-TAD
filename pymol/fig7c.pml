# fig7c.pml -- Fig 7C (and optional 7B): naked duplexes through one cohesin ring, full-atom model
# usage (full model loaded):   @fig7c.pml
#   then render:   view fig7C_side      png fig7c_side.png, width=6559, height=3531, ray=1
#                  view fig7C_through   png fig7c_through.png, width=6559, height=3531, ray=1
#   optional wide view of the pinch (set WIDE = 1 below, rerun):
#                  view fig7B_pinch     png fig7b_pinch.png, width=6559, height=3531, ray=1
#
# After the script has run, the ring can be changed interactively (each takes a second or two):
#   ring_size 55            set the ring radius in A (centre of the tube)
#   ring_move 5, 0, -10     shift the ring in A: sideways (u), sideways (v), along the ring axis
#   ring_tilt 10, 0         tilt the ring in degrees about u and v
#   ring_fit                re-fit position and tilt at the current radius
#   ring_fit min            find the smallest radius that threads only the two naked duplexes
#   ring_info               print radius, clearance, and every DNA strand passing through the ring
# The script also prints a measurement table for all three tethers.
# Restore the rest of the model afterwards with:  delete cohesin_ring   and   show_as lines, all
python
import json, itertools, numpy as np
from pymol import cmd
from pymol.cgo import CYLINDER, SPHERE, COLOR
TETHER = "T_I5"      # tether to draw: T_I5, T_I7 or T_I2
REP = "lines"        # representation for the shown atoms; use the one you used for Figs 6 and 7C
WIDE = 0             # 1 = also show +-8 kb of both arms and store the view fig7B_pinch
RING_R = 60.0        # starting ring radius in A (centre of the tube); a duplex is 20 A wide
RING_TUBE = 10.0     # tube radius in A
FIT_MIN = 1          # 1 = start from the smallest radius that threads only the two naked duplexes
MEASURE_ALL = 1      # 1 = print the measurement table for all three tethers (slower)
THROUGH_SLAB = 60.0  # half-thickness (A) of the slice shown in the through-the-ring view; larger shows more depth
DATA = {"T_I5": {"centre": [2605.0726352766233, -1807.6708677971678, -114.33986548278418], "axis": [0.22422871049121929, 0.9740756316460084, -0.029969137873397823], "R": 200.0, "tube": 12.0, "linkers": [{"kind": "cohesin", "naked": [13277028, 13277450], "span": [13276494, 13277984], "nucs": [2594, 2599], "wide_span": [13269239, 13285239], "wide_nucs": [2556, 2637]}, {"kind": "mirror", "naked": [12924189, 12924935], "span": [12923677, 12925481], "nucs": [746, 751], "wide_span": [12916562, 12932562], "wide_nucs": [711, 787]}]}, "T_I7": {"centre": [-1517.2190098856647, -2032.5996715608835, -1875.2857102520368], "axis": [-0.5907727191553939, -0.8025296123012136, 0.08326953633475949], "R": 200.0, "tube": 12.0, "linkers": [{"kind": "cohesin", "naked": [13368429, 13369041], "span": [13367895, 13369575], "nucs": [3070, 3075], "wide_span": [13360735, 13376735], "wide_nucs": [3032, 3113]}, {"kind": "mirror", "naked": [13005675, 13006287], "span": [13005141, 13006821], "nucs": [1169, 1174], "wide_span": [12997981, 13013981], "wide_nucs": [1131, 1212]}]}, "T_I2": {"centre": [2003.2747859177505, -3849.8819118361844, -526.2328559825509], "axis": [0.8218303041783047, 0.19856809505228126, 0.5340090474528529], "R": 100.0, "tube": 12.0, "linkers": [{"kind": "cohesin", "naked": [12783706, 12784370], "span": [12783188, 12784872], "nucs": [53, 58], "wide_span": [12776038, 12792038], "wide_nucs": [21, 93]}, {"kind": "mirror", "naked": [13354368, 13354980], "span": [13353834, 13355514], "nucs": [2999, 3004], "wide_span": [13346674, 13362674], "wide_nucs": [2961, 3042]}]}}
A = lambda g: g - 12771345
B = lambda g: 13423346 - g
def dna(a, b): return "((chain A and resi %d-%d) or (chain B and resi %d-%d))" % (A(a), A(b), B(b), B(a))
def hist(i0, i1): return "(chain " + "+".join("P%d*" % i for i in range(i0, i1 + 1)) + ")"
def unit(x): x = np.asarray(x, float); return x / np.linalg.norm(x)
def xyz(sel): return cmd.get_coords(sel) if cmd.count_atoms(sel) else np.zeros((0, 3))
def mind(s1, s2, cut):
    a = xyz("(%s) within %g of (%s)" % (s1, cut, s2)); b = xyz("(%s) within %g of (%s)" % (s2, cut, s1))
    if len(a) == 0 or len(b) == 0: return None
    return float(np.sqrt(((a[:, None] - b[None]) ** 2).sum(-1)).min()), cmd.count_atoms("(%s) within 3 of (%s)" % (s1, s2))
def axis_pts(a0, b0):
    m, mb = {}, {}
    cmd.iterate_state(1, "chain A and resi %d-%d and name C1'" % (A(a0), A(b0)), "m[int(resi)] = (x, y, z)", space={"m": m})
    cmd.iterate_state(1, "chain B and resi %d-%d and name C1'" % (B(b0), B(a0)), "mb[652001 - int(resi)] = (x, y, z)", space={"mb": mb})
    P = np.array([(np.array(m[r]) + np.array(mb[r])) / 2 for r in sorted(m) if r in mb])
    if len(P) < 12: return np.zeros((0, 3))
    return np.array([P[i:i + 10].mean(0) for i in range(len(P) - 10)])
def ring_frame(ax):
    u = unit(np.cross(ax, [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0])); return u, np.cross(ax, u)
def torus_gap(P, c, ax, R, tube):
    q = P - c; h = q @ ax; rho = np.linalg.norm(q - np.outer(h, ax), axis=1)
    return np.sqrt((rho - R) ** 2 + h ** 2) - tube
def crossings(S, c, ax):
    out = []
    if len(S) < 2: return out
    h = (S - c) @ ax
    for i in np.where(np.sign(h[:-1]) != np.sign(h[1:]))[0]:
        q = S[i] + h[i] / (h[i] - h[i + 1]) * (S[i + 1] - S[i]); out.append((np.linalg.norm((q - c) - np.dot(q - c, ax) * ax), q))
    return out
class Ring:
    """state of one tether's ring: frame, the two naked duplex axes, nearby atoms"""
    def __init__(self, t):
        self.c0 = np.array(t["centre"], float); self.a0 = unit(t["axis"]); self.u, self.v = ring_frame(self.a0)
        self.S = [axis_pts(*l["naked"]) for l in t["linkers"]]
        cmd.pseudoatom("_rc", pos=[float(x) for x in self.c0])
        self.P = xyz("(chain A+B or chain P*) within 260 of _rc")
        self.D = xyz("(chain A+B) within 260 of _rc")
        naked = " or ".join(dna(*l["naked"]) for l in t["linkers"])
        self.Q = xyz("((chain A+B or chain P*) within 260 of _rc) and not (%s)" % naked); cmd.delete("_rc")
        self.p = [0, 0, 0, 0, 0]; self.R = RING_R
    def place(self, p=None):
        p = self.p if p is None else p
        ax = unit(self.a0 + np.tan(np.radians(p[0])) * self.u + np.tan(np.radians(p[1])) * self.v)
        return self.c0 + p[2] * self.u + p[3] * self.v + p[4] * self.a0, ax
    def check(self, p=None, R=None):
        """returns (ok, clearance, crossings per duplex, number of other atoms inside the opening)"""
        R = self.R if R is None else R; c, ax = self.place(p)
        inner = R - RING_TUBE
        cr = [crossings(S, c, ax) for S in self.S]
        through = [[z for z in x if z[0] < inner] for x in cr]
        graze = [[z for z in x if inner <= z[0] < R + RING_TUBE + 10] for x in cr]
        ok = all(len(x) == 1 and x[0][0] + 10.0 + 2.0 <= inner for x in through) and not any(graze)
        # any other DNA or protein (not the two naked duplexes) inside the opening
        q = self.Q - c; h = q @ ax; rho = np.linalg.norm(q - np.outer(h, ax), axis=1)
        others = int(((np.abs(h) < 8.0) & (rho < inner)).sum())
        extra = sum(max(len(x) - 1, 0) for x in through) + sum(1 for x in through if len(x) == 0) + sum(len(x) for x in graze)
        self.last_bad = extra * 1000 + others
        ok = ok and others == 0
        clr = float(torus_gap(self.P, c, ax, R, RING_TUBE).min()) if len(self.P) else 99.0
        return ok, clr, through, others
    def fit(self, R=None, quick=False):
        R = self.R if R is None else R; best = None
        grid = ((-20, 0, 20), (-20, 0, 20), (-30, 0, 30), (-30, 0, 30), (-30, 0, 30)) if quick else \
               ((-20, -10, 0, 10, 20), (-20, -10, 0, 10, 20), (-40, -20, 0, 20, 40), (-40, -20, 0, 20, 40), (-40, -20, 0, 20, 40))
        fallback = None
        for cand in itertools.product(*grid):
            ok, clr, _, _ = self.check(cand, R)
            if ok and (best is None or clr > best[0]): best = (clr, list(cand))
            key = (self.last_bad, -clr)
            if fallback is None or key < fallback[0]: fallback = (key, list(cand))
        if best is None:
            self.R, self.p = R, fallback[1]; return False
        for step in ((5, 5, 10, 10, 10), (2.5, 2.5, 5, 5, 5), (1, 1, 2, 2, 2)):
            for d in itertools.product((-1, 0, 1), repeat=5):
                cand = [b + k * s_ for b, k, s_ in zip(best[1], d, step)]
                ok, clr, _, _ = self.check(cand, R)
                if ok and clr > best[0]: best = (clr, cand)
        self.R, self.p = R, best[1]; return True
    def fit_min(self):
        for R in (40, 45, 50, 55, 60, 70, 80, 90, 100, 120):
            if self.fit(R, quick=True):
                self.fit(R); return True
        self.fit(RING_R); return False
def draw(ring):
    c, ax = ring.place(); u, v = ring_frame(ax); col = (0.20, 0.20, 0.19); o = []
    P = [c + ring.R * (np.cos(x) * u + np.sin(x) * v) for x in np.linspace(0, 2 * np.pi, 145)]
    for p, q in zip(P[:-1], P[1:]): o += [CYLINDER, *p, *q, RING_TUBE, *col, *col, COLOR, *col, SPHERE, *q, RING_TUBE]
    cmd.delete("cohesin_ring"); cmd.load_cgo(o, "cohesin_ring")
def info(ring, name=""):
    ok, clr, through, others = ring.check()
    rr = ["%.0f" % x[0][0] if len(x) == 1 else "%d crossings" % len(x) for x in through]
    print("%s ring radius %.0f A (opening %.0f A across) | clearance %.1f A | cohesin-site duplex: %s A, mirror duplex: %s A from centre | other atoms in the opening: %d | %s"
          % (name, ring.R, 2 * (ring.R - RING_TUBE), clr, rr[0], rr[1], others, "OK" if ok else "CHECK: not only the two naked duplexes pass through"))
cmd.set("suspend_updates", 1)
RINGS = {}
# ---------- measurements for all three tethers ----------
fits = []
for tid, t in DATA.items():
    if not MEASURE_ALL and tid != TETHER: continue
    ring = Ring(t); RINGS[tid] = ring
    if any(len(S) == 0 for S in ring.S): fits.append((tid, None)); continue
    good = ring.fit_min() if FIT_MIN else ring.fit()
    fits.append((tid, good))
print("\n%-5s %6s %9s %14s %14s %10s %10s" % ("ring", "R (A)", "axis sep", "DNA-DNA", "DNA-hist", "hist-ring", "DNA-ring"))
for tid, good in fits:
    t = DATA[tid]; ring = RINGS[tid]
    if good is None: print("%-5s tether atoms not loaded" % tid); continue
    c, ax = ring.place(); ok, clr, through, others = ring.check()
    sep = np.linalg.norm(through[0][0][1] - through[1][0][1]) if all(len(x) == 1 for x in through) else float("nan")
    nk = [dna(*l["naked"]) for l in t["linkers"]]; naked = "(%s or %s)" % (nk[0], nk[1])
    dd = mind(nk[0], nk[1], 12); dh = mind(naked, "chain P*", 12)
    hr = float(torus_gap(xyz("chain P* within 200 of (%s)" % naked), c, ax, ring.R, RING_TUBE).min()) if cmd.count_atoms("chain P* within 200 of (%s)" % naked) else None
    dr = float(torus_gap(ring.D, c, ax, ring.R, RING_TUBE).min()) if len(ring.D) else None
    f = lambda x: ">12" if x is None else "%.1f (%d)" % x
    g = lambda x: ">40" if x is None else "%.1f" % x
    print("%-5s %6.0f %9.1f %14s %14s %10s %10s" % (tid, ring.R, sep, f(dd), f(dh), g(hr), g(dr)))
print("axis sep = distance between the two duplex axes where they cross the ring plane. Other columns = closest atom distance in A,")
print("in brackets the number of atoms closer than 3 A. DNA-hist = naked duplexes against any histone atom.")
print("hist-ring / DNA-ring = closest atom to the surface of the ring tube (negative means inside the tube)")
for tid, good in fits:
    if good is None: continue
    if not good: print("%s: no position found that threads only the two naked duplexes; showing the best attempt, check it" % tid)
    info(RINGS[tid], tid)
print("")
# ---------- draw the chosen tether ----------
t = DATA[TETHER]; ring = RINGS.get(TETHER) or Ring(t); RINGS[TETHER] = ring
parts = []
for l in t["linkers"]:
    nm = "site_" + l["kind"]
    cmd.select(nm + "_dna", dna(*l["naked"]))
    span, nucs = (l["wide_span"], l["wide_nucs"]) if WIDE else (l["span"], l["nucs"])
    cmd.select(nm + "_ctx", dna(*span) + " or " + hist(*nucs))
    parts.append(nm)
cmd.select("tether_all", " or ".join(p + "_ctx" for p in parts))
cmd.hide("everything", "not tether_all")
cmd.show_as(REP, "tether_all")
cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("depth_cue", 0)
cmd.set("orthoscopic", 1); cmd.set("field_of_view", 20)
cmd.set_color("t_hist", [0.84, 0.82, 0.78]); cmd.set_color("t_dna", [0.62, 0.60, 0.56])
cmd.set_color("t_coh", [0.122, 0.435, 0.471])     # cohesin-site duplex: teal
cmd.set_color("t_mir", [0.851, 0.451, 0.122])     # mirror (inverted-copy) duplex: orange
cmd.color("t_dna", "tether_all and chain A+B"); cmd.color("t_hist", "tether_all and chain P*")
cmd.color("t_coh", "site_cohesin_dna"); cmd.color("t_mir", "site_mirror_dna")
draw(ring)
# ---------- views ----------
def store_views():
    c, ax = ring.place(); ok, clr, through, _ = ring.check()
    if all(len(x) == 1 for x in through): s = through[1][0][1] - through[0][0][1]
    else: s = ring_frame(ax)[0]
    s = unit(s - np.dot(s, ax) * ax)
    def setv(right, up, d, name, depth):
        z = np.cross(right, up)
        cmd.set_view((right[0], up[0], z[0], right[1], up[1], z[1], right[2], up[2], z[2], 0.0, 0.0, -d, c[0], c[1], c[2], d - depth, d + depth, -20.0))
        cmd.view(name, "store")
        print("%-14s camera distance %6.0f A  ->  %.3f px per A at 3531 px height" % (name, d, 3531 / (2 * d * np.tan(np.radians(10)))))
    dist = max(500.0, ring.R * 6.0)
    setv(s, np.cross(ax, s), dist, "fig7C_through", THROUGH_SLAB)
    setv(ax, s, dist * 1.4, "fig7C_side", 900.0)
    if WIDE:
        S0 = ring.S[0]; arm = unit(S0[-1] - S0[0]); arm = unit(arm - np.dot(arm, s) * s)
        setv(arm, s, 6000.0, "fig7B_pinch", 4000.0)
store_views()
cmd.set("suspend_updates", 0)
print("stored views: fig7C_side, fig7C_through" + (", fig7B_pinch" if WIDE else ""))
# ---------- interactive ring commands ----------
def _redraw():
    draw(ring); info(ring, TETHER); store_views()
def ring_size(r):
    ring.R = float(r); _redraw()
def ring_move(du=0, dv=0, da=0):
    ring.p[2] += float(du); ring.p[3] += float(dv); ring.p[4] += float(da); _redraw()
def ring_tilt(a=0, b=0):
    ring.p[0] += float(a); ring.p[1] += float(b); _redraw()
def ring_fit(mode=""):
    good = ring.fit_min() if str(mode).strip() == "min" else ring.fit()
    if not good: print("no position found that threads only the two naked duplexes at this radius; try a larger ring_size")
    _redraw()
def ring_info():
    info(ring, TETHER)
for f_ in (ring_size, ring_move, ring_tilt, ring_fit, ring_info): cmd.extend(f_.__name__, f_)
print("ring commands ready: ring_size R | ring_move du, dv, da | ring_tilt a, b | ring_fit [min] | ring_info")
python end
