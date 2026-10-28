# fig5_style_views.pml -- Fig 5 on the loaded final model: colours AND the three Fig 5 views, all in this one file
# Run after the model is loaded:   @fig5_style_views.pml
# then  view fig5A / view fig5B / view fig5C  and ray each one yourself.
# Fog: depth cue on with PyMOL's default fog_start (0.45); adjust with e.g.  set fog_start, 0.6
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
set field_of_view, 20
set suspend_updates, off
# the three concentric views (same rotation and centre, only the camera distance differs)
set_view (-0.122136, 0.832715, 0.540059, -0.508779, 0.414676, -0.754440, -0.852188, -0.366915, 0.373025, 0.0, 0.0, -892.85, 231.497941, 549.481079, -421.883972, 257.03, 1526.12, -20.0)
view fig5C, store
set_view (-0.122136, 0.832715, 0.540059, -0.508779, 0.414676, -0.754440, -0.852188, -0.366915, 0.373025, 0.0, 0.0, -3203.37, 231.497941, 549.481079, -421.883972, -1263.57, 7667.75, -20.0)
view fig5B, store
set_view (-0.122136, 0.832715, 0.540059, -0.508779, 0.414676, -0.754440, -0.852188, -0.366915, 0.373025, 0.0, 0.0, -22904.29, 231.497941, 549.481079, -421.883972, 16622.47, 29183.55, -20.0)
view fig5A, store
view fig5A
print("Fig 5 ready: view fig5A / view fig5B / view fig5C, then ray")
