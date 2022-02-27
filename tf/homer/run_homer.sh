#!/usr/bin/sh


 reuse -q .homer-4.10

#Make input bed files for emVars and for all MPRA variants, with 200 bp window
awk '$6=="TGWAS" && $58=="Enhancer_Skew"' /dir/T_cell_MPRA/mpra_annot/mpra_data_merge.txt | awk '{print $1"\t"($57-100)"\t"($2+100)}' > MPRA.emVAR_200bp_pad.bed

awk '$6=="TGWAS"' /dir/T_cell_MPRA/mpra_annot/mpra_data_merge.txt  | awk '{print $1"\t"($57-100)"\t"($2+100)}' > MPRA.all_200bp_pad.bed

#Run HOMER enrichment, using all MPRA variants as the background and emVars as the foreground
w=200 #window parameter
findMotifsGenome.pl /dir/T_cell_MPRA/mpra_annot/MPRA.emVAR_${w}bp_pad.bed hg19 mpra_emvar_${w} -size ${w} -bg /dir/T_cell_MPRA/mpra_annot/MPRA.all_${w}bp_pad.bed
