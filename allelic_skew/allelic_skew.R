library(readxl)
library(data.table)
library(stringr)
library(IRanges)
library(SummarizedExperiment)

dir<-"/dir/T_cell_MPRA/"
##########################################################
###Plot ATAC-seq allelic skew against MPRA allelic skew###

#Read in allelic skew data from Calderon et al. (PMID: 31570894), Supplementary Table 1, "significant_ASCs" tab
dat.skew<-data.frame(read_excel("/Users/michaelguo/Dropbox/JR_MPRA/asc/41588_2019_505_MOESM3_ESM.xlsx", sheet=10), stringsAsFactors = F)
dat.skew<-dat.skew[!grepl("\\-S", dat.skew$cell_type),] #remove stimulated cell types

#Parse ATAC allelic skew file
dat.skew$chr<-str_split_fixed(dat.skew$het_id, "\\_", 5)[,2]
dat.skew$pos<-as.numeric(str_split_fixed(dat.skew$het_id, "\\_", 5)[,3])
dat.skew$chr_pos<-paste(dat.skew$chr, dat.skew$pos, sep=":")

#Collapse allelic skew counts across cell types for each SNP
dat.skew.collapse<-data.frame(chr_pos=unique(dat.skew$chr_pos), het_id=0,atac_ref=0, atac_alt=0, stringsAsFactors = F)
for(i in 1:nrow(dat.skew.collapse)){
  dat.skew.collapse[i,]$atac_ref<-sum(subset(dat.skew, chr_pos==dat.skew.collapse[i,]$chr_pos)$refCount)
  dat.skew.collapse[i,]$atac_alt<-sum(subset(dat.skew, chr_pos==dat.skew.collapse[i,]$chr_pos)$altCount)
  dat.skew.collapse[i,]$het_id<-subset(dat.skew, chr_pos==dat.skew.collapse[i,]$chr_pos)$het_id[1]
}

#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra,project=="TGWAS", select=c(chr, snp_start, snp_end, LogSkew, mpra_sig))

#Merge ATAC allelic skew data with MPRA data
mpra$chr_pos<-paste0(mpra$chr, ":", mpra$snp_end)
mpra<-merge(mpra, dat.skew.collapse, by="chr_pos", all.x=T, all.y=F)
mpra$allelic_skew<-mpra$atac_alt/(mpra$atac_alt+mpra$atac_ref)

#Plot MPRA Skew against ATAC allelic skew
mpra.skew<-subset(mpra, mpra_sig%in%c("Enhancer_Skew" ,"Enhancer_nSkew") & !is.na(het_id) )
plot(x=mpra.skew$allelic_skew, y=mpra.skew$LogSkew, pch=16, cex=1.2 ,col="darkgray", xlim=c(0,1), ylim=c(-1,1))
points(x=mpra.skew[mpra.skew$mpra_sig=="Enhancer_Skew",]$allelic_skew, y=mpra.skew[mpra.skew$mpra_sig=="Enhancer_Skew",]$LogSkew, pch=16, cex=1.2 ,col="red", xlim=c(0,1), ylim=c(-1,1))
abline(v=0.5, lty=2, col="blue")
abline(h=0, lty=2, col="blue")
abline(lm(mpra.skew$LogSkew~mpra.skew$allelic_skew))

#Generate correlation statistics
summary(lm(mpra.skew$LogSkew~mpra.skew$allelic_skew))


#Plot by MPRA activity
dat.plot<-data.frame( mpra_sig=c("nEnhancer_nSkew", "Enhancer_nSkew", "Enhancer_Skew"), proportion=0,stringsAsFactors = F)
for(i in 1:nrow(dat.plot)){
  dat.plot[i,]$proportion<-nrow(subset(mpra, !is.na(allelic_skew) & mpra_sig==dat.plot[i,]$mpra_sig))/nrow(subset(mpra, mpra_sig==dat.plot[i,]$mpra_sig))
}
dat.plot$mpra_sig<-factor(dat.plot$mpra_sig, levels=c("nEnhancer_nSkew", "Enhancer_nSkew", "Enhancer_Skew"))

#Plot proportion of pCREs and emVars that show allelic skew in Calderon ATAC-seq data
ggplot(dat.plot, aes(fill=mpra_sig, y=proportion, x=mpra_sig)) + 
  geom_bar(position="dodge", stat="identity", color="black")+
  theme_bw()+scale_fill_brewer(palette = "Blues")+coord_cartesian(ylim=c(0,0.02))+
  theme(axis.text.x = element_text( angle = 45, hjust=1))+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")

#Calculate enrichment p-values of MPRA variants in allelic skew variants
#For pCREs
prop.test(x = c(nrow(subset(mpra, !is.na(allelic_skew) & mpra_sig=="Enhancer_nSkew")), nrow(subset(mpra, !is.na(allelic_skew) & mpra_sig=="nEnhancer_nSkew"))), 
          n =c(nrow(subset(mpra,  mpra_sig=="Enhancer_nSkew")), nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew"))))
#For emVars
prop.test(x = c(nrow(subset(mpra, !is.na(allelic_skew) & mpra_sig=="Enhancer_Skew")), nrow(subset(mpra, !is.na(allelic_skew) & mpra_sig=="nEnhancer_nSkew"))), 
          n =c(nrow(subset(mpra,  mpra_sig=="Enhancer_Skew")), nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew"))))

