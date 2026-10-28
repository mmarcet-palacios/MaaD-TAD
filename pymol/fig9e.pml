# fig9e.pml -- the D3/D4 arms side by side in the full-atom model (chr1:13.18-13.40 Mb)
# Start a FRESH PyMOL (nothing loaded), then:
#     cd MaaD-TAD/pymol
#     @fig9e.pml
# Loads blocks 42-64 one at a time and keeps only chr1:13,184,371-13,395,000 (nucleosomes 2109-3209),
# so memory stays low (~4-5 M atoms kept). Takes a few minutes. Then colours the arms, draws dashed
# rungs every 5 kb between matching positions (145.75 kb apart) and stores view "d3d4".
#
# Colours
#   context (loop between the arms, histones)  light grey
#   arm A, D3+D4  13,184,371-13,224,371        sky blue
#   arm C, D3+D4  13,330,121-13,370,121        deep blue
#   Micro-C line  13,225,000-13,250,000 and 13,370,000-13,395,000   orange (both arms)
#   rungs                                      black dashes (CGO, so they ray-trace)
# Residue numbering: chain A resi = genomic - 12,771,345, chain B resi = 13,423,346 - genomic

reinitialize
set suspend_updates, on
python
from pymol import cmd
import os
BLOCKS = r"../blocks_final"
OFF_A, OFF_B = 12771345, 13423346
G0, G1 = 13184371, 13395000                  # kept window
N0, N1 = 2109, 3209                          # nucleosomes whose dyads fall in the window
keep_dna = "(chain A and resi %d-%d) or (chain B and resi %d-%d)" % (G0 - OFF_A, G1 - OFF_A, OFF_B - G1, OFF_B - G0)
for b in range(N0 // 50, N1 // 50 + 1):
    f = os.path.join(BLOCKS, "block_%02d.cif" % b)
    name = "b%02d" % b
    cmd.load(f, name)
    lo, hi = max(N0, 50 * b), min(N1, 50 * b + 49)
    hist = " or ".join("chain P%d*" % i for i in range(lo, hi + 1)) if lo <= hi else "none"
    cmd.remove("%s and not ((%s) or (%s))" % (name, keep_dna, hist))
    print("block %02d: kept %d atoms" % (b, cmd.count_atoms(name)))
cmd.create("d3d4", "b*")
cmd.delete("b*")
python end

bg_color white
set ray_opaque_background, 1
set depth_cue, 0
set_color ctx_dna, [0.78, 0.77, 0.74]
set_color ctx_hist, [0.88, 0.87, 0.84]
set_color armA, [0.53, 0.78, 0.92]
set_color armC, [0.12, 0.33, 0.62]
set_color line_or, [0.85, 0.45, 0.12]
hide everything
show_as lines, d3d4
color ctx_hist, d3d4 and chain P*
color ctx_dna, d3d4 and chain A+B
# arm A (D3+D4) and arm C (D3+D4)
color armA, d3d4 and ((chain A and resi 413026-453026) or (chain B and resi 198975-238975))
color armC, d3d4 and ((chain A and resi 558776-598776) or (chain B and resi 53225-93225))
# the 25 kb seen in Micro-C, on both arms
color line_or, d3d4 and ((chain A and resi 453655-478655) or (chain B and resi 173346-198346))
color line_or, d3d4 and ((chain A and resi 598655-623655) or (chain B and resi 28346-53346))

# dashed rungs every 5 kb between matching positions on the two arms (offset 145,750 bp)
python
from pymol import cmd
from pymol.cgo import CYLINDER
import numpy as np
OFF_A, OFF_B, SHIFT = 12771345, 13423346, 145750
def centre(g, half=250):
    sel = "d3d4 and ((chain A and resi %d-%d) or (chain B and resi %d-%d))" % (
        g - half - OFF_A, g + half - OFF_A, OFF_B - g - half, OFF_B - g + half)
    xyz = cmd.get_coords(sel)
    return None if xyz is None else xyz.mean(0)
obj, R, DASH, GAP = [], 8.0, 35.0, 30.0
for g in range(13185000, 13250001, 5000):
    p, q = centre(g), centre(g + SHIFT)
    if p is None or q is None:
        continue
    v = q - p; L = float(np.linalg.norm(v)); u = v / L; s = 0.0
    while s < L:
        a = p + u * s; e = p + u * min(s + DASH, L)
        obj += [CYLINDER, *map(float, a), *map(float, e), R, 0, 0, 0, 0, 0, 0]
        s += DASH + GAP
    print("rung %.3f Mb <-> %.3f Mb: %.0f A" % (g / 1e6, (g + SHIFT) / 1e6, L))
cmd.load_cgo(obj, "rungs")
python end

orient d3d4 and chain A+B
set suspend_updates, off
view d3d4, store
# Rotate until the two arms lie side by side, then:   view d3d4, store
# final render (same width as Fig 5):
#   view d3d4
#   png ../figures/final/fig9e.png, width=6559, height=4920, ray=1
