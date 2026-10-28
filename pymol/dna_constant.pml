# dna_constant.pml -- one neutral colour for DNA, paler histones, so the tubes carry the colour-coding
set suspend_updates, on
set_color dna_ink, [0.420, 0.404, 0.361]
set_color histone_pale, [0.851, 0.831, 0.784]
color dna_ink, chain A+B
color histone_pale, chain P*
set suspend_updates, off
