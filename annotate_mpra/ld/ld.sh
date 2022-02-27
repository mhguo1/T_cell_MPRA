#!/usr/bin/sh

#This script identifies all LD proxies for MPRA variants
cut -f2 /dir/T_cell_MPRA/tag_analysis/CMS_OL12_JURKAT_emVAR_20200827.out  | tail -n+2 > mpra.snp.list

for chr in {1..22} X
do
plink2 \
--vcf /ref/1kg/EUR.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.recode.vcf.gz \
--r2 \
--ld-snp-list mpra.snp.list \
--ld-window-kb 2000 \
--ld-window 999999 \
--ld-window-r2 0.2 \
--out mpra.chr${chr}

rm -f *nosex
rm -f *log
done

#Note that 1000 Genomes Phase 3 vcfs were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ and subsetted to 503 individuals of European ancestry. The ID column in the vcf was then recoded to CHROM:POS:REF:ALT
