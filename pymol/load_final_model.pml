# load_final_model.pml -- load the FINAL model (all 68 blocks of blocks_final) for Figs 4A, 5 and 7C
# Fresh PyMOL:   cd MaaD-TAD/pymol
#                @load_final_model.pml        (about 5-20 min, the window stays frozen while it loads)
# Then style one figure at a time:
#   @fig4_stage5.pml                      Fig 4A stage 5 (grey, camera of stages 2-4)
#   @fig5_style_views.pml                 Fig 5 (blue DNA, orange histones), then  view fig5A / fig5B / fig5C
#   @fig7c.pml  then  ring_size 45        Fig 7C (T_I5), then  view fig7C_side
# Ray each one yourself, e.g.   png ../figures/final/fig5A_ray.png, width=6559, height=3531, ray=1
reinitialize
set auto_dss, 0
set suspend_updates, on
python
from pymol import cmd
import os, time
t = time.time()
for i in range(68):
    cmd.load(os.path.join(r"../blocks_final", "block_%02d.cif" % i), "block_%02d" % i)
cmd.show_as("lines", "all")
cmd.set("suspend_updates", 0)
print("final model loaded: 68 blocks, %d atoms, %.0f s" % (cmd.count_atoms("all"), time.time() - t))
python end
