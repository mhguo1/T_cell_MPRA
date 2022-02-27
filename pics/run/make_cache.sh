#!/usr/bin/sh

#Make cache file for PICS
cat /dir/T_cell_MPRA/annotate_mpra/ld/*.ld |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | awk 'BEGIN {OFS=FS="\t"} {if (NR>1) {print "chr"$3, "chr"$6, "-", $7, $3":"$6}}' | awk 'BEGIN {OFS=FS="\t"} {gsub(/:/, "\t", $5) }1' | awk 'BEGIN {OFS=FS="\t"} { if ($1!="CHR_A") {print $1, $2, $3, $4, $7","$11}}' | awk '$1!=$2' | grep -v SNP | awk '$4>0.8' | python flip_cache.py | python remove_cache_dups.py > cache

#Make cache index file for PICS
cat /dir/T_cell_MPRA/annotate_mpra/ld/*ld | grep -v ^SNP |  awk '{print "chr"$3}' | sort | uniq > cacheindex

#Create an emtpy file for missing SNPs for PICS
touch cachemissing

