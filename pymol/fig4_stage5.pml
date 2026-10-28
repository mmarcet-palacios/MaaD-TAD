# fig4_stage5.pml -- Fig 4A stage 5: the full-atom model in grey, same orientation as stages 2-4
# 1. Fresh PyMOL, load all 68 blocks (the whole locus is in this panel). No other scripts (not fig7c.pml).
# 2. @fig4_stage5.pml        (sets representation, colours and camera; renders nothing)
# 3. Quick look, seconds:    png fig4_stage5_test.png, width=600, height=600, ray=0
# 4. Final:                  png fig4_stage5.png, width=1800, height=1800, ray=1
#    Do not touch the mouse between steps 2 and 4: the camera must stay where this script put it.
set suspend_updates, on
show_as lines, all
# colours: dark DNA, pale histones (same as Fig 4B and dna_constant.pml)
set_color dna_ink, [0.420, 0.404, 0.361]
set_color histone_pale, [0.851, 0.831, 0.784]
color dna_ink, chain A+B
color histone_pale, chain P*
bg_color white
set ray_opaque_background, 1
set depth_cue, 1
# camera: the projection used for stages 2-4 (centre = mean nucleosome position, 40,000 A away)
set_view (\
    -0.529923,    0.789860,   -0.308712,\
    -0.606986,   -0.607489,   -0.512372,\
    -0.592241,   -0.084134,    0.801356,\
     0.000000,    0.000000, -40000.000000,\
   -80.395282,   29.261712,  101.881969,\
 28000.000000, 52000.000000,  -20.000000 )
set orthoscopic, on
set field_of_view, 20
set suspend_updates, off
