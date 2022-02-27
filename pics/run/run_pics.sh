#!/usr/bin/sh
ld=0.8
for prefix in Celiac Crohns IBD PBC RA SLE T1D UC MS Psoriasis
do
while read line
do
snp=`echo $line | awk '{print $1}'` #snp should start with chr:pos:ref:alt and should start with "chr"
pval=`echo $line | awk '{print $2}'`

perl pics.pl $snp $pval /dir/T_cell_MPRA/annotate_mpra/ld | grep "^chr"

done < /dir/T_cell_MPRA/gwas/gwas_catalog/${prefix}.catalog_sumstats.txt > /dir/T_cell_MPRA/pics/run/${prefix}.ld_${ld}.pics 2> /dir/T_cell_MPRA/pics/run/${prefix}.ld_${ld}.pics.err
done
