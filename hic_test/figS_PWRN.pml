# figS_PWRN.pml -- supplementary figure, gap-free test model PWRN_chr15 (panels B and C)
# 1. Fresh PyMOL:  cd MaaD-TAD/hic_test   then   @figS_PWRN.pml
#    (loads 33 blocks, a few minutes; renders nothing)
# 2. Panel B:  view figS_B   quick look:  png figS_B_test.png, width=600, height=600, ray=0
#    final:   png figS_PWRN_B.png, width=3000, height=3000, ray=1
#    scale (same for both loci): at 3000 px, 100 nm = 162.1 px; 500 nm = 810 px
#    labels are the object labels_B (disable labels_B to render without them)
# 3. Panel C:  figS_closeup   (cartoon of the close-up nucleosomes only, then view figS_C)
#    final:   png figS_PWRN_C.png, width=3000, height=3000, ray=1
#    Do not touch the mouse between a view command and its png.
reinitialize
set auto_dss, 0
set suspend_updates, on
python
from pymol import cmd
import os, time
t = time.time()
for i in range(33):
    cmd.load(os.path.join(r"PWRN_chr15_H1/blocks_final", "block_%02d.cif" % i), "block_%02d" % i)
print("PWRN_chr15: %d atoms, %.0f s" % (cmd.count_atoms("all"), time.time() - t))
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
# L4 arm 1 24,196,001-24,236,000: 211 nucleosomes
color armA, (chain A and resi 66001-106000) or (chain B and resi 204001-244000)
color armA, chain P347*+P348*+P349*+P350*+P351*+P352*+P353*+P354*+P355*+P356*+P357*+P358*+P359*+P360*+P361*+P362*+P363*+P364*+P365*+P366*+P367*+P368*+P369*+P370*+P371*+P372*+P373*+P374*+P375*+P376*+P377*+P378*+P379*+P380*+P381*+P382*+P383*+P384*+P385*+P386*+P387*+P388*+P389*+P390*+P391*+P392*+P393*+P394*+P395*+P396*+P397*+P398*+P399*+P400*+P401*+P402*+P403*+P404*+P405*+P406*+P407*+P408*+P409*+P410*+P411*+P412*+P413*+P414*+P415*+P416*+P417*+P418*+P419*+P420*+P421*+P422*+P423*+P424*+P425*+P426*+P427*+P428*+P429*+P430*+P431*+P432*+P433*+P434*+P435*+P436*+P437*+P438*+P439*+P440*+P441*+P442*+P443*+P444*+P445*+P446*+P447*+P448*+P449*+P450*+P451*+P452*+P453*+P454*+P455*+P456*+P457*+P458*+P459*+P460*+P461*+P462*+P463*+P464*+P465*+P466*+P467*+P468*+P469*+P470*+P471*+P472*+P473*+P474*+P475*+P476*+P477*+P478*+P479*+P480*+P481*+P482*+P483*+P484*+P485*+P486*+P487*+P488*+P489*+P490*+P491*+P492*+P493*+P494*+P495*+P496*+P497*+P498*+P499*+P500*+P501*+P502*+P503*+P504*+P505*+P506*+P507*+P508*+P509*+P510*+P511*+P512*+P513*+P514*+P515*+P516*+P517*+P518*+P519*+P520*+P521*+P522*+P523*+P524*+P525*+P526*+P527*+P528*+P529*+P530*+P531*+P532*+P533*+P534*+P535*+P536*+P537*+P538*+P539*+P540*+P541*+P542*+P543*+P544*+P545*+P546*+P547*+P548*+P549*+P550*+P551*+P552*+P553*+P554*+P555*+P556*+P557*
# L4 arm 2 24,294,001-24,342,000: 252 nucleosomes
color armC, (chain A and resi 164001-212000) or (chain B and resi 98001-146000)
color armC, chain P863*+P864*+P865*+P866*+P867*+P868*+P869*+P870*+P871*+P872*+P873*+P874*+P875*+P876*+P877*+P878*+P879*+P880*+P881*+P882*+P883*+P884*+P885*+P886*+P887*+P888*+P889*+P890*+P891*+P892*+P893*+P894*+P895*+P896*+P897*+P898*+P899*+P900*+P901*+P902*+P903*+P904*+P905*+P906*+P907*+P908*+P909*+P910*+P911*+P912*+P913*+P914*+P915*+P916*+P917*+P918*+P919*+P920*+P921*+P922*+P923*+P924*+P925*+P926*+P927*+P928*+P929*+P930*+P931*+P932*+P933*+P934*+P935*+P936*+P937*+P938*+P939*+P940*+P941*+P942*+P943*+P944*+P945*+P946*+P947*+P948*+P949*+P950*+P951*+P952*+P953*+P954*+P955*+P956*+P957*+P958*+P959*+P960*+P961*+P962*+P963*+P964*+P965*+P966*+P967*+P968*+P969*+P970*+P971*+P972*+P973*+P974*+P975*+P976*+P977*+P978*+P979*+P980*+P981*+P982*+P983*+P984*+P985*+P986*+P987*+P988*+P989*+P990*+P991*+P992*+P993*+P994*+P995*+P996*+P997*+P998*+P999*+P1000*+P1001*+P1002*+P1003*+P1004*+P1005*+P1006*+P1007*+P1008*+P1009*+P1010*+P1011*+P1012*+P1013*+P1014*+P1015*+P1016*+P1017*+P1018*+P1019*+P1020*+P1021*+P1022*+P1023*+P1024*+P1025*+P1026*+P1027*+P1028*+P1029*+P1030*+P1031*+P1032*+P1033*+P1034*+P1035*+P1036*+P1037*+P1038*+P1039*+P1040*+P1041*+P1042*+P1043*+P1044*+P1045*+P1046*+P1047*+P1048*+P1049*+P1050*+P1051*+P1052*+P1053*+P1054*+P1055*+P1056*+P1057*+P1058*+P1059*+P1060*+P1061*+P1062*+P1063*+P1064*+P1065*+P1066*+P1067*+P1068*+P1069*+P1070*+P1071*+P1072*+P1073*+P1074*+P1075*+P1076*+P1077*+P1078*+P1079*+P1080*+P1081*+P1082*+P1083*+P1084*+P1085*+P1086*+P1087*+P1088*+P1089*+P1090*+P1091*+P1092*+P1093*+P1094*+P1095*+P1096*+P1097*+P1098*+P1099*+P1100*+P1101*+P1102*+P1103*+P1104*+P1105*+P1106*+P1107*+P1108*+P1109*+P1110*+P1111*+P1112*+P1113*+P1114*
pseudoatom labels_B, pos=[-427.3, 4134.2, 1638.5], label="L4 arm 1"
pseudoatom labels_B, pos=[-3619.3, -1113.9, 1773.2], label="L4 arm 2"
set label_size, -450
set label_color, black
set label_font_id, 7
hide nonbonded, labels_B
set_view (-0.042351, 0.999028, 0.012204, 0.919621, 0.034205, 0.391314, 0.390516, 0.027795, -0.920176, 0.000000, 0.000000, -52483.410467, -476.090109, -55.812682, 7.401891, 43828.563693, 61138.257240, -20.000000)
view figS_B, store
# close-up: nucleosomes 448-455 (dyads 24,215,236-24,216,566) in L4 arm 1
select closeup, (chain A and resi 85163-86639) or (chain B and resi 223362-224838) or chain P448*+P449*+P450*+P451*+P452*+P453*+P454*+P455*
set_view (0.515410, 0.821089, 0.245286, -0.095428, -0.229460, 0.968629, 0.851614, -0.522649, -0.039911, 0.000000, 0.000000, -2602.803939, -2231.626356, 4625.060393, 2494.909770, 2346.207112, 2859.400765, -20.000000)
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
