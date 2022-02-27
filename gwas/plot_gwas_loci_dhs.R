library(stringr)
library(IRanges)
library(SummarizedExperiment)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(readxl)

#############################
dir<-"/dir/T_cell_MPRA"

gwas.order<- c("Crohns","IBD", "MS","Psoriasis","RA","T1D","UC" )

#Read in file with allele directions.
#This file contains one row for each emVar
#"risk_allele" column shows whether the ref (A) or alt (B) allele is the risk allele
#Remaining columns show for each disease, whether the alt (B) allele increases or decreases risk for that disease. If not associated, then it's NA
dat.alleles<-data.frame(read_excel(paste0(dir, "/gwas/mpra_allele_annotation.xlsx"), sheet=1), stringsAsFactors = F)
dat.alleles<-subset(dat.alleles, select=c(label,risk_allele, MS_risk_allele, Psoriasis_risk_allele,RA_risk_allele,T1D_risk_allele,UC_risk_allele,IBD_risk_allele,Crohns_risk_allele))
names(dat.alleles)<-gsub("_risk_allele", "", names(dat.alleles))

###Add in PICS fine-mapping data###

#Read in MPRA data and subset for emVars in DHS
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
#mpra<-read.delim("/Users/michaelguo/Dropbox/JR_MPRA/tables/mpra_data_merge.txt", header=T, stringsAsFactors = F, sep="\t")
#names(mpra)<-gsub("_CS_", "_PP_", names(mpra))
mpra$label<-paste0(mpra$SNP, "_", mpra$rsid)
mpra$dhs<-mpra$dhs_Tcell_merged

mpra<-mpra[mpra$project=="TGWAS" & mpra$dhs==1 & mpra$mpra_sig=="Enhancer_Skew",]

#For each disease, add in a new column. 
#If a given MPRA variant is in the 95% credible set, then assign as 1. 
#If it's a proxy SNP (r2>0.8), but not in credible set, then arbitrarily assign to 0.5
mpra$Crohns<-0
mpra$Crohns[mpra$Crohns_pics>=0.95 | (mpra$Crohns_pics>=0.01 & mpra$Crohns_PP_running<=0.95)]<-1
mpra$Crohns[!is.na(mpra$Crohns_pics)&mpra$Crohns_pics<0.95 & (mpra$Crohns_pics<0.01| mpra$Crohns_PP_running>0.95)]<-0.5
mpra$IBD<-0
mpra$IBD[mpra$IBD_pics>=0.95 | (mpra$IBD_pics>=0.01 & mpra$IBD_PP_running<=0.95)]<-1
mpra$IBD[!is.na(mpra$IBD_pics)&mpra$IBD_pics<0.95 & (mpra$IBD_pics<0.01| mpra$IBD_PP_running>0.95)]<-0.5
mpra$MS<-0
mpra$MS[mpra$MS_pics>=0.95 | (mpra$MS_pics>=0.01 & mpra$MS_PP_running<=0.95)]<-1
mpra$MS[!is.na(mpra$MS_pics)&mpra$MS_pics<0.95 & (mpra$MS_pics<0.01| mpra$MS_PP_running>0.95)]<-0.5
mpra$Psoriasis<-0
mpra$Psoriasis[mpra$Psoriasis_pics>=0.95 | mpra$Psoriasis_pics>=0.01 & mpra$Psoriasis_PP_running<=0.95]<-1
mpra$Psoriasis[!is.na(mpra$Psoriasis_pics)& mpra$Psoriasis_pics<0.95 &(mpra$Psoriasis_pics<0.01| mpra$Psoriasis_PP_running>0.95)]<-0.5
mpra$RA<-0
mpra$RA[mpra$RA_pics>=0.95 | (mpra$RA_pics>=0.01 & mpra$RA_PP_running<=0.95)]<-1
mpra$RA[!is.na(mpra$RA_pics)&mpra$RA_pics<0.95 & (mpra$RA_pics<0.01| mpra$RA_PP_running>0.95)]<-0.5
mpra$T1D<-0
mpra$T1D[mpra$T1D_pics>=0.95 | (mpra$T1D_pics>=0.01 & mpra$T1D_PP_running<=0.95)]<-1
mpra$T1D[!is.na(mpra$T1D_pics)&mpra$T1D_pics<0.95 & (mpra$T1D_pics<0.01| mpra$T1D_PP_running>0.95)]<-0.5
mpra$UC<-0
mpra$UC[mpra$UC_pics>=0.95 | (mpra$UC_pics>=0.01 & mpra$UC_PP_running<=0.95)]<-1
mpra$UC[!is.na(mpra$UC_pics)& mpra$UC_pics<0.95 &(mpra$UC_pics<0.01| mpra$UC_PP_running>0.95)]<-0.5
mpra<-subset(mpra, select=c(SNP, label, LogSkew, Skew.logFDR, MS, Psoriasis, RA, T1D, UC, IBD, Crohns))


mpra.reshape<-melt(data=mpra, id.vars=c("SNP", "label", "LogSkew", "Skew.logFDR"), measure.vars = names(mpra)[5:ncol(mpra)])
mpra.reshape$FC<-2^(abs(mpra.reshape$LogSkew))


#Plot allelic fold change for emVars
mpra.reshape<-mpra.reshape[order(mpra.reshape$label),]
mpra.skew<-mpra.reshape[!duplicated(mpra.reshape$label),]
mpra.skew<-mpra.skew[order(-mpra.skew$FC),]
snp.order<-unique(mpra.skew$label)
mpra.skew$label <-factor(mpra.skew$label, levels=snp.order)

ggplot(data=mpra.skew, aes(x=label, y=FC,label=label))+
  geom_bar(stat="identity",aes(fill=Skew.logFDR), color="black")+scale_fill_gradient(limits=c(1,1.4),low = "grey",high = "#A80707")+theme_bw()+
  theme(axis.text.x=element_text(color = "black", angle=45, hjust=1, size=8),axis.text.y=element_text(color = "black", size=8))



#Determine directionality of GWAS effect for each emVar with respect to the MPRA activity-increasing allele
#NB: Will need to manually fix rs654690 in Illustrator
mpra.skew<-merge(mpra.skew, dat.alleles[,c("label", "risk_allele")], by="label", all.x=T, all.y=T)
mpra.skew$direction<-0
for(i in 1:nrow(mpra.skew)){
  if(mpra.skew[i,]$risk_allele=="A"){
    if(mpra.skew[i,]$LogSkew>0){
      mpra.skew[i,]$direction<- -1
    }else{
      mpra.skew[i,]$direction<- 1
    }
  }else{
    if(mpra.skew[i,]$LogSkew>0){
      mpra.skew[i,]$direction<- 1
    }else{
      mpra.skew[i,]$direction<- -1
    }
  }
}
ggplot(data=mpra.skew, aes(x=label,y=0,fill=direction, label=SNP))+
  geom_tile(color="black")+theme_bw()+theme(legend.position="n")+scale_fill_gradient(low ="royalblue2" ,high ="lightcoral")+
  theme(axis.text.x=element_text(color = "black",angle=45, hjust=1, size=8),axis.text.y=element_text(color = "black", size=8))



#Plot membership of each emVar in credible sets for each disease
#Plot in dark purple if in PICS credible set
#Plot in light purple if it is in LD with a GWAS hit (r2>0.8), but not in credible set.
dat.pics<-mpra.reshape
dat.pics$gwas<-dat.pics$variable
dat.pics$pics<-dat.pics$value

dat.pics$gwas<-factor(dat.pics$gwas, levels=rev(gwas.order))
dat.pics$label<-factor(dat.pics$label, levels=snp.order)

ggplot(data=dat.pics, aes(x=label, y=gwas, label=SNP))+theme_bw()+
  geom_tile(color="black",aes(fill=pics))+scale_fill_gradient(low = "white",high = "purple")+
  theme(axis.text.x=element_text(color = "black",angle=45, hjust=1, size=8),axis.text.y=element_text(color = "black", size=8))
