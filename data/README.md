# data/ — pipeline inputs (CC BY 4.0)

| File | Contents | Used by |
|---|---|---|
| 1KX5.cif | nucleosome core particle, 147 bp on a histone octamer (PDB 1KX5) | Builder |
| b2_onelayer.pse | one layer of the two-start zigzag built from 1KX5 (PyMOL session) | analyse_b2.py, via build_zigzag.py (Placer, Builder) |
| gap-filled.xyz | CSynth path of the locus, 326 beads of 2 kb (bead start, x, y, z), gap filled | Folder (refold) |
| gap_sections.txt | sections of the 50 kb assembly gap | Folder, Placer |
| nucleosomes_GRCh38_resolved.bed | 3,376 nucleosome positions (490 measured, GM12878 MNase-seq lifted to GRCh38; the rest placeholders at 190 bp) | Folder, Placer, Fig 1A |
| nucleosomes_tethered.bed | the same after opening naked DNA at the three tethers (3,359 nucleosomes) | Placer |
| homology_restraints.json | the 16 homology blocks (arms, orientation, identity) | Folder, Fig 1C |
| ctcf_cohesin_mirror_sites.csv | CTCF/cohesin sites and their predicted mirror sites (Table 1; reference table) | |
| guide_refolded.xyz, guide_tethered.xyz | Folder outputs (refolded path; path with pairing and tethers) | Placer |
| tethers.json, tethers_built.json | cohesin ring geometry and waypoints (current; as built) | Builder, Checker |
| verify_nucleosome_source.py | lifts the 490 measured dyads back to hg19 and checks them against the UCSC track (needs pyliftover, internet) | |

Locus: GRCh38 chr1:12,771,346–13,423,345 (sequence in ../locus_GRCh38.fa).
