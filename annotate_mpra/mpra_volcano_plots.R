library(stringr)
library(ggplot2)


dir<-"/dir/T_cell_MPRA"
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra, project=="TGWAS" & !is.na(LogSkew),select=c(SNP, A.log2FC, B.log2FC, LogSkew,mpra_sig))


#For each MPRA variant, find the fold change of the allele with the largest absolute fold change
mpra$max_abs_fc<-NA
for(i in 1:nrow(mpra)){
  if(abs(mpra[i,]$A.log2FC)>abs(mpra[i,]$B.log2FC)){
    mpra[i,]$max_abs_fc<-mpra[i,]$A.log2FC
  }else{
    mpra[i,]$max_abs_fc<-mpra[i,]$B.log2FC
  }
}

#Create temporary dataframe with all MPRA alleles
mpra.all<-mpra
mpra.all$mpra_sig_expanded<-"all"

#Create temporary dataframe with only pCRE alleles
mpra.pcre<-subset(mpra, mpra_sig%in%c("Enhancer_nSkew","Enhancer_Skew"))
mpra.pcre$mpra_sig_expanded<-"pcre"

#Create temporary dataframe with only emVars
mpra.emvar<-subset(mpra, mpra_sig=="Enhancer_Skew")
mpra.emvar$mpra_sig_expanded<-"emvar"

#Combined temporary dataframes
mpra.expanded<-rbind(mpra.all, mpra.pcre, mpra.emvar)


#Generate volcano plot 
ggplot() + 
  geom_point(data=mpra,aes(x=max_abs_fc, y=LogSkew), alpha=0.5, pch=16, color="grey50", size=1) +
  geom_point(data=mpra[mpra$mpra_sig!="nEnhancer_nSkew",],aes(x=max_abs_fc, y=LogSkew), alpha=0.2, pch=16, color="goldenrod", size=0.5) +
  geom_point(data=mpra[mpra$mpra_sig=="Enhancer_Skew",],aes(x=max_abs_fc, y=LogSkew), alpha=1, pch=16, color="firebrick4", size=1) +
  theme_bw()+ylab("MPRA Log2(allelic skew)")+xlab("Expression log2(max_allele)")+
  coord_cartesian(xlim=c(-2,7), ylim=c(-2,2))


#Generate density plot for MPRA signal (absolute foldchange) 
mpra.expanded$mpra_sig_expanded <- factor(mpra.expanded$mpra_sig_expanded, levels=c("all", "pcre", "emvar"))
ggplot(mpra.expanded, aes(max_abs_fc, fill=mpra_sig_expanded))+geom_density(alpha=0.4, col="white")+
  scale_fill_manual(values=c("grey50","goldenrod", "firebrick4"))+theme_bw()+coord_cartesian(xlim=c(-2,7))

#Generate density plot for MPRA allelic skew
ggplot(mpra.expanded, aes(LogSkew, fill=mpra_sig_expanded))+geom_density(alpha=0.4, col="white")+
  scale_fill_manual(values=c("grey50","goldenrod", "firebrick4"))+theme_bw()+coord_cartesian(xlim=c(-2,2))
