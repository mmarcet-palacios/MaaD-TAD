# figS_IGHC.pml -- supplementary figure, gap-free test model IGHC_chr14 (panels B and C)
# 1. Fresh PyMOL:  cd MaaD-TAD/hic_test   then   @figS_IGHC.pml
#    (loads 23 blocks, a few minutes; renders nothing)
# 2. Panel B:  view figS_B   quick look:  png figS_B_test.png, width=600, height=600, ray=0
#    final:   png figS_IGHC_B.png, width=3000, height=3000, ray=1
#    scale (same for both loci): at 3000 px, 100 nm = 162.1 px; 500 nm = 810 px
#    labels are the object labels_B (disable labels_B to render without them)
# 3. Panel C:  figS_closeup   (cartoon of the close-up nucleosomes only, then view figS_C)
#    final:   png figS_IGHC_C.png, width=3000, height=3000, ray=1
#    Do not touch the mouse between a view command and its png.
reinitialize
set auto_dss, 0
set suspend_updates, on
python
from pymol import cmd
import os, time
t = time.time()
for i in range(23):
    cmd.load(os.path.join(r"IGHC_chr14_H1/blocks_final", "block_%02d.cif" % i), "block_%02d" % i)
print("IGHC_chr14: %d atoms, %.0f s" % (cmd.count_atoms("all"), time.time() - t))
python end
show_as lines, all
bg_color white
set ray_opaque_background, 1
set depth_cue, 1
set orthoscopic, on
set field_of_view, 20
set_color hp_dna, [0.620, 0.600, 0.560]
set_color hp_hist, [0.840, 0.820, 0.780]
set_color g_gene, [0.184, 0.561, 0.357]
set_color g_pseudo, [0.300, 0.300, 0.280]
set_color armA, [0.530, 0.780, 0.920]
set_color armC, [0.120, 0.330, 0.620]
color hp_dna, chain A+B
color hp_hist, chain P*
# IGHA2 105,583,731-105,588,395: 25 nucleosomes
color g_gene, (chain A and resi 43731-48395) or (chain B and resi 161606-166270)
color g_gene, chain P230*+P231*+P232*+P233*+P234*+P235*+P236*+P237*+P238*+P239*+P240*+P241*+P242*+P243*+P244*+P245*+P246*+P247*+P248*+P249*+P250*+P251*+P252*+P253*+P254*
# IGHE 105,597,691-105,601,728: 21 nucleosomes
color g_gene, (chain A and resi 57691-61728) or (chain B and resi 148273-152310)
color g_gene, chain P304*+P305*+P306*+P307*+P308*+P309*+P310*+P311*+P312*+P313*+P314*+P315*+P316*+P317*+P318*+P319*+P320*+P321*+P322*+P323*+P324*
# IGHG4 105,620,506-105,626,066: 29 nucleosomes
color g_gene, (chain A and resi 80506-86066) or (chain B and resi 123935-129495)
color g_gene, chain P424*+P425*+P426*+P427*+P428*+P429*+P430*+P431*+P432*+P433*+P434*+P435*+P436*+P437*+P438*+P439*+P440*+P441*+P442*+P443*+P444*+P445*+P446*+P447*+P448*+P449*+P450*+P451*+P452*
# IGHG2 105,639,559-105,644,790: 27 nucleosomes
color g_gene, (chain A and resi 99559-104790) or (chain B and resi 105211-110442)
color g_gene, chain P524*+P525*+P526*+P527*+P528*+P529*+P530*+P531*+P532*+P533*+P534*+P535*+P536*+P537*+P538*+P539*+P540*+P541*+P542*+P543*+P544*+P545*+P546*+P547*+P548*+P549*+P550*
# IGHGP 105,664,633-105,669,843: 27 nucleosomes
color g_pseudo, (chain A and resi 124633-129843) or (chain B and resi 80158-85368)
color g_pseudo, chain P656*+P657*+P658*+P659*+P660*+P661*+P662*+P663*+P664*+P665*+P666*+P667*+P668*+P669*+P670*+P671*+P672*+P673*+P674*+P675*+P676*+P677*+P678*+P679*+P680*+P681*+P682*
# IGHA1 105,703,995-105,708,665: 24 nucleosomes
color g_gene, (chain A and resi 163995-168665) or (chain B and resi 41336-46006)
color g_gene, chain P863*+P864*+P865*+P866*+P867*+P868*+P869*+P870*+P871*+P872*+P873*+P874*+P875*+P876*+P877*+P878*+P879*+P880*+P881*+P882*+P883*+P884*+P885*+P886*
# IGHEP1 105,721,794-105,722,527: 4 nucleosomes
color g_pseudo, (chain A and resi 181794-182527) or (chain B and resi 27474-28207)
color g_pseudo, chain P956*+P957*+P958*+P959*
# IGHG1 105,736,343-105,743,071: 35 nucleosomes
color g_gene, (chain A and resi 196343-203071) or (chain B and resi 6930-13658)
color g_gene, chain P1033*+P1034*+P1035*+P1036*+P1037*+P1038*+P1039*+P1040*+P1041*+P1042*+P1043*+P1044*+P1045*+P1046*+P1047*+P1048*+P1049*+P1050*+P1051*+P1052*+P1053*+P1054*+P1055*+P1056*+P1057*+P1058*+P1059*+P1060*+P1061*+P1062*+P1063*+P1064*+P1065*+P1066*+P1067*
pseudoatom labels_B, pos=[2414.0, 1720.4, 314.4], label="IGHA2"
pseudoatom labels_B, pos=[-188.3, 2348.0, 550.8], label="IGHE"
pseudoatom labels_B, pos=[-3989.8, -1118.7, -718.6], label="IGHG4"
pseudoatom labels_B, pos=[-5334.0, -1002.0, 726.0], label="IGHG2"
pseudoatom labels_B, pos=[-1594.0, -3600.4, 1154.5], label="IGHGP"
pseudoatom labels_B, pos=[2578.2, 74.8, -1263.0], label="IGHA1"
pseudoatom labels_B, pos=[-859.3, 1983.5, -1302.7], label="IGHEP1"
pseudoatom labels_B, pos=[-2879.0, -2259.9, -1748.8], label="IGHG1"
set label_size, -450
set label_color, black
set label_font_id, 7
hide nonbonded, labels_B
set_view (-0.800024, -0.560629, -0.213674, 0.176466, -0.560267, 0.809296, -0.573430, 0.609750, 0.547159, 0.000000, 0.000000, -52483.410467, -31.420713, -291.786971, 671.268430, 45285.838720, 59680.982214, -20.000000)
view figS_B, store
# close-up: nucleosomes 1046-1053 (dyads 105,738,902-105,740,232) in IGHG1
select closeup, (chain A and resi 198829-200305) or (chain B and resi 9696-11172) or chain P1046*+P1047*+P1048*+P1049*+P1050*+P1051*+P1052*+P1053*
set_view (-0.824357, -0.520389, 0.222781, -0.518742, 0.852005, 0.070675, -0.226589, -0.057304, -0.972303, 0.000000, 0.000000, -2650.435042, -2306.455788, -1805.013786, -2156.149097, 2394.191239, 2906.678845, -20.000000)
view figS_C, store
python
def figS_closeup():
    cmd.hide("everything"); cmd.disable("labels_B")
    cmd.set("cartoon_ring_mode", 3); cmd.set("cartoon_ring_finder", 2)
    cmd.show("cartoon", "closeup"); cmd.view("figS_C", "recall")
cmd.extend("figS_closeup", figS_closeup)
python end
deselect
view figS_B
set suspend_updates, off
