# render_all.pml -- Figs 4A (stage 5), 5A-C and 6C, redrawn with thicker lines and lighter fog, from the FINAL model (blocks_final)
# Run in the PyMOL window (Schroedinger PyMOL, the GUI, not headless), from this folder:
#     cd MaaD-TAD/pymol
#     @render_all.pml
# Large views use the graphics card (draw), which takes seconds per image, and ray-tracing 52 M atoms ran out of memory.
# Fig 7C is made separately (fig7c_render.pml, open-source PyMOL, no licence watermark). Outputs go to ../figures/final/ (nothing else is overwritten).
reinitialize
set max_threads, 24
set hash_max, 255
# faster loading of 52 M atoms: no secondary-structure assignment, no screen updates until all blocks are in
set auto_dss, 0
set suspend_updates, on
python
from pymol import cmd
import os, time
BLOCKS = r"../blocks_final"
t = time.time()
for i in range(68):
    cmd.load(os.path.join(BLOCKS, "block_%02d.cif" % i), "block_%02d" % i)
    print("  block %02d loaded, %.0f s" % (i, time.time() - t), flush=True)
cmd.set("suspend_updates", 0)
print("loaded 68 blocks, %d atoms, %.0f s" % (cmd.count_atoms("all"), time.time() - t), flush=True)
def stamp(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)
cmd.extend("stamp", stamp)
python end

# drawn (OpenGL) lines are 1 px wide at 6559 px and the fog is stronger than in ray-tracing:
# thicker lines and a later fog start bring the look back to the ray-traced panels
set line_width, 2.5
set fog_start, 0.65
set line_smooth, 1

# ---- Fig 4A stage 5 (grey, camera of stages 2-4) ----
stamp Fig4 stage5 start
@fig4_stage5.pml
draw 1800, 1800, antialias=2
png ../figures/final/fig4_stage5_v2.png
stamp Fig4 stage5 done

# ---- Fig 5A-C (DNA blue, histones orange, lines, orthographic, depth cue) ----
stamp Fig5 start
set suspend_updates, on
show_as lines, all
set_color f6_dna, [0.247, 0.435, 0.659]
set_color f6_hist, [0.851, 0.451, 0.122]
color f6_dna, chain A+B
color f6_hist, chain P*
bg_color white
set ray_opaque_background, 1
set depth_cue, 1
set orthoscopic, on
set suspend_updates, off
@fig5_views.pml
view fig5A
draw 6559, 3531, antialias=2
png ../figures/final/fig5A_v2.png
stamp Fig5A done
view fig5B
draw 6559, 3531, antialias=2
png ../figures/final/fig5B_v2.png
stamp Fig5B done
view fig5C
draw 6559, 3531, antialias=2
png ../figures/final/fig5C_v2.png
stamp Fig5C done

# ---- Fig 6C (I11/I12 hairpin, mirror-image genes) ----
stamp Fig6C start
show_as lines, all
@fig6c.pml
view fig6C
draw 6559, 3531, antialias=2
png ../figures/final/fig6c_v2.png
stamp Fig6C done
stamp model figures done

print("ALL DONE")
