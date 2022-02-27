library(stringr)

gwas.order<- c("Crohns","MS","Psoriasis", "RA","T1D","UC", "IBD" )

dir<-"/dir/T_cell_MPRA"

###Histogram of number of MPRA hits per locus
#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra, project=="TGWAS",select=c( chr, snp_end,ld_snp, mpra_sig))

#Read in and merge LD files across chromosomes
dat.ld.all<-data.frame(lead_snp=character(), ld_snp=character(), r2=numeric(), stringsAsFactors = F)
for(c in c(1:22, "X")){
  dat.ld.chr<-read.delim(paste(dir, "/annotate_mpra/ld/mpra.chr", c, ".ld", sep=""), header=T, stringsAsFactors = F, sep="")
  dat.ld.chr<-subset(dat.ld.chr, select=c(SNP_A, SNP_B, R2))
  names(dat.ld.chr)<-c("lead_snp","ld_snp", "r2")
  dat.ld.all<-rbind(dat.ld.all, dat.ld.chr)
  rm(dat.ld.chr)
}
mpra.merge<-merge(dat.ld.all, mpra, by="ld_snp", all.x=F, all.y=F)

#Remove MHC region and also 17:44073889:A:G as they have very high number of MPRA variants tested
dat.sig<-subset(mpra.merge, !(chr==6 & snp_end>29691116 & snp_end<33054976) & !(lead_snp%in%c("17:44073889:A:G")) & r2>=0.8 & mpra_sig=="Enhancer_Skew" )
dat.plot<-data.frame(table(dat.sig$lead_snp), stringsAsFactors = F)

hist(table(dat.sig$lead_snp), xlim=c(0,8),breaks=seq(0,8,1), ylim=c(0,120), col="skyblue2", main="Histogram of # of MPRA SNPs per locus", xlab="# of MPRA SNPs", ylab="# of Loci")



####Plot Number of loci explained,###
#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra, project=="TGWAS")

#Create dataframe with one row per disease. Then, count total number of loci tested, as well as # of loci with at least one emVar
dat.disease<-data.frame(disease=gwas.order, total=0, mpra=0, stringsAsFactors = F)
for(i in 1:nrow(dat.disease)){
  d<-dat.disease[i,]$disease
  temp<-mpra[!is.na(mpra[,which(names(mpra)==paste0(d, "_pval"))]) & mpra[,which(names(mpra)==paste0(d, "_pval"))]>=-log10(5e-8),]
  dat.disease[i,]$total<-length(unique(temp$lead_snp))
  dat.disease[i,]$mpra<-length(unique(temp[temp$mpra_sig=="Enhancer_Skew",]$lead_snp))
}

dat.disease<-dat.disease[order(-dat.disease$total),]
dat.disease$total<-dat.disease$total-dat.disease$mpra
dat.plot<-melt(dat.disease, id.vars="disease")
dat.plot$disease<-factor(dat.plot$disease, levels=unique(dat.disease$disease))

ggplot(data=dat.plot[dat.plot$variable!="mpra_atac",], aes(x=disease, y=value,  fill =variable)) +geom_bar(stat="identity")+
  theme_bw()+scale_fill_brewer(palette = "Set2")+ylab("# of loci")+
  theme(axis.text.x=element_text(color = "black",angle=45,hjust=1))
