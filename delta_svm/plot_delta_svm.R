library(data.table)
library(stringr)

dir<-"/dir/T_cell_MPRA"

#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra,project=="TGWAS" & mpra_sig=="Enhancer_Skew", select=c(chr, snp_start, snp_end, SNP, mpra_sig, LogSkew))

dat.deltasvm.cd4<-read.delim(paste0(dir, "/delta_svm/mpra_snps_E2_naive_CD4_deltaSVM.txt"), header=F, stringsAsFactors = F, sep="\t")
names(dat.deltasvm.cd4)<-c("SNP", "delta_svm")
mpra.merge<-merge(mpra, dat.deltasvm.cd4, by="SNP", all.x=F, all.y=F)

#Plot MPRA allele logSkew against deltaSVM score
ggplot(mpra.merge, aes(x=delta_svm, y=LogSkew))+geom_point(color="red", fill="red",pch=16, alpha=0.5, size=2)+
  theme_bw()+ 
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+ 
  geom_abline(intercept = as.numeric(lm(mpra.merge$LogSkew~mpra.merge$delta_svm)$coefficients[1]), slope = as.numeric(lm(mpra.merge$LogSkew~mpra.merge$delta_svm)$coefficients[2]), color="blue")+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-2.5, 2.5))

#Generate correlation statistics
summary(lm(mpra.merge$delta_svm~mpra.merge$LogSkew))
