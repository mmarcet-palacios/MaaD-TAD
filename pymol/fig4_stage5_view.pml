# fig4_stage5_view.pml -- the Fig 4A view (same orientation as stages 2-4) for the full-atom model
# Run after the model is loaded:  @fig4_stage5_view.pml
set orthoscopic, on
set_view (-0.529923, 0.789860, -0.308712, -0.606986, -0.607489, -0.512372, -0.592241, -0.084134, 0.801356, 0.000000, 0.000000, -40000.000000, -80.395282, 29.261712, 101.881969, 28000.000000, 52000.000000, 1.000000)
bg_color white
# neutral colours as in the strip (comment out to keep your current colours)
@dna_constant.pml
# square image for the strip; ray-tracing 55 M atoms can take a while -- use ray=0 for a quick draw
png fig4_stage5.png, width=1800, height=1800, dpi=300, ray=1
