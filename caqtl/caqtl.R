library(readxl)
library(data.table)
library(stringr)
library(IRanges)
library(SummarizedExperiment)
library(data.table)

dir<-"/dir/T_cell_MPRA"

##########################################################
###Plot caQTL skew against MPRA allelic skew##############

#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra,project=="TGWAS", select=c(chr, snp_start, snp_end, LogSkew, mpra_sig))

#Read in Gate et al ATAC-QTL data (PMID: 29988122), Supplemental Table 6
dat.qtl<-data.frame(read_excel("/Users/michaelguo/Dropbox/JR_MPRA/qtl/41588_2018_156_MOESM8_ESM.xlsx", sheet=2), stringsAsFactors = F)
dat.qtl$chr<-str_split_fixed(dat.qtl$SNP, "\\:", 2)[,1]
dat.qtl$snp_end<-str_split_fixed(str_split_fixed(dat.qtl$SNP, "\\_", 2)[,1], "\\:", 2)[,2]
dat.qtl<-subset(dat.qtl, select=c(chr, snp_end, beta, p.value))
names(dat.qtl)[3:4]<-c("atac_qtl_beta", "atac_qtl_pval")

#Merge in caQTL data with MPRA data
mpra<-merge(mpra, dat.qtl, by=c("chr", "snp_end"), all.x=T, all.y=F)

#Plot proportion of MPRA variants that are caQTLs
dat.plot<-data.frame(table(subset(mpra, !is.na(atac_qtl_pval))$mpra_sig)/table(mpra$mpra_sig),stringsAsFactors = F)
names(dat.plot)<-c("mpra_sig", "proportion")
dat.plot$mpra_sig<-factor(dat.plot$mpra_sig, levels=c("nEnhancer_nSkew", "Enhancer_nSkew", "Enhancer_Skew"))
dat.plot$mpra_sig<-factor(dat.plot$mpra_sig, levels=c("nEnhancer_nSkew", "Enhancer_nSkew", "Enhancer_Skew"))

ggplot(dat.plot, aes(fill=mpra_sig, x=mpra_sig, y=proportion)) + 
  geom_bar(position="dodge", stat="identity", color="black")+
  theme_bw()+scale_fill_brewer(palette = "Blues")+coord_cartesian(ylim=c(0,0.02))+
  theme(axis.text.x = element_text( angle = 45, hjust=1))+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")


#Calculate enrichment p-values of MPRA variants in allelic skew variants
#For pCREs
prop.test(x = c(nrow(subset(mpra, !is.na(atac_qtl_pval) & mpra_sig=="Enhancer_nSkew")), nrow(subset(mpra, !is.na(atac_qtl_pval)& mpra_sig=="nEnhancer_nSkew"))), 
          n =c(nrow(subset(mpra,  mpra_sig=="Enhancer_nSkew")), nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew"))))
#For emVars
prop.test(x = c(nrow(subset(mpra, !is.na(atac_qtl_pval) & mpra_sig=="Enhancer_Skew")), nrow(subset(mpra, !is.na(atac_qtl_pval)& mpra_sig=="nEnhancer_nSkew"))), 
          n =c(nrow(subset(mpra,  mpra_sig=="Enhancer_Skew")), nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew"))))


#Plot MPRA Skew against caQTL skew 
mpra.qtl<-subset(mpra, !is.na(atac_qtl_pval) & mpra_sig%in%c("Enhancer_Skew", "Enhancer_nSkew"))
plot(x=mpra.qtl$atac_qtl_beta, y=mpra.qtl$LogSkew, pch=16, cex=1.2 ,col="darkgray", xlim=c(0.3,0.7), ylim=c(-1.5,1.5))
points(x=mpra.qtl[mpra.qtl$mpra_sig=="Enhancer_Skew",]$atac_qtl_beta, y=mpra.qtl[mpra.qtl$mpra_sig=="Enhancer_Skew",]$LogSkew, pch=16, cex=1.2 ,col="red", xlim=c(0,1), ylim=c(-1,1))
abline(v=0.5, lty=2, col="blue")
abline(h=0, lty=2, col="blue")
abline(lm(mpra.qtl$LogSkew~mpra.qtl$atac_qtl_beta))






