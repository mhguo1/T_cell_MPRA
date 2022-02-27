###T1D analysis--Onegut-Gumuscu et al Fine-mapping analyses
library(readxl)
library(stringr)

dir<-"/dir/T_cell_MPRA"

#Read in MPRA data
dat.mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")

#Subset for loci with at least one emVar 
emvar_loci_lead_snps<-unique(subset(dat.mpra, mpra_sig=="Enhancer_Skew")$lead_snp) #find lead SNPs of loci w at least one emVar
dat.mpra<-subset(dat.mpra, lead_snp%in%emvar_loci_lead_snps) #Subset MPRA data to only loci with at least one emVar
emvar_loci_lead_rsid<-unique(subset(dat.mpra, lead_snp==ld_snp)$rsid) #Generate a list of rsIDs for these loci

#Read in finemap data
#Download Table S1 from Onegut-Gumuscu et al: https://static-content.springer.com/esm/art%3A10.1038%2Fng.3245/MediaObjects/41588_2015_BFng3245_MOESM391_ESM.xls
dat.finemap<-data.frame(read_excel("41588_2015_BFng3245_MOESM391_ESM.xls", sheet=1, skip=1), stringsAsFactors = F)
dat.finemap<-subset(dat.finemap, !is.na(name.rs))

#Assign posterior probabilities for proxy SNPs to their corresponding credible set SNPs
for(i in 1:nrow(dat.finemap)){
  if(dat.finemap[i,]$snp.type=="proxy"){
    dat.finemap[i,]$pp<-max(dat.finemap[dat.finemap$name.rs==dat.finemap[i,]$cred.snp.rs,]$pp)
  }
}
dat.finemap<-subset(dat.finemap, select=c(name.rs, index.snp.rs,snp.type, pp, r2))

#Find lead SNP IDs and rsIDs for T1D loci that are also tested in the MPRA as lead SNPs
lead_t1d_rsid<-unique(dat.finemap[dat.finemap$index.snp.rs%in%emvar_loci_lead_rsid,]$index.snp.rs)
lead_t1d_snpid<-unique(subset(dat.mpra, rsid%in%lead_t1d_rsid)$lead_snp)

#Test enrichment of MPRA SNPs in T1D credible sets
dat.finemap$pp<-as.numeric(dat.finemap$pp)
dat.enrichment<-data.frame(pp=c(0.01, 0.05, 0.1,0.15,0.2, 0.25,0.3,0.4,0.5 ),  
                     a=0, b=0, c=0, d=0,fold=0, p=0, stringsAsFactors = F)

for(i in 1:nrow(dat.enrichment)){
  dat.finemap.subset<-subset(dat.finemap, pp>=dat.enrichment[i,]$pp &index.snp.rs%in%lead_t1d_rsid) #Generate list of fine-mapped SNPs
  dat.mpra<-subset(dat.mpra, lead_snp%in%lead_t1d_snpid) #Extract these from the MPRA data
  a<-nrow(subset(dat.mpra, mpra_sig=="Enhancer_Skew" & rsid%in%dat.finemap.subset$name.rs)) #Find # of emVars that are fine-mapped
  b<-nrow(subset(dat.mpra, mpra_sig=="Enhancer_Skew" & !(rsid%in%dat.finemap.subset$name.rs))) #Find # of emVars that are not fine-mapped
  c<-nrow(subset(dat.mpra, mpra_sig!="Enhancer_Skew" & rsid%in%dat.finemap.subset$name.rs)) #Find # of non-emVars that are fine-mapped
  d<-nrow(subset(dat.mpra, mpra_sig!="Enhancer_Skew" & !(rsid%in%dat.finemap.subset$name.rs))) #Find # of non-emVars that are not fine-mapped
  dat.enrichment[i,]$a<-a
  dat.enrichment[i,]$b<-b
  dat.enrichment[i,]$c<-c
  dat.enrichment[i,]$d<-d
  dat.enrichment[i,]$fold<-(a/(a+b))/(c/(c+d)) #Calculate enrichment
  dat.enrichment[i,]$p<-fisher.test(rbind(c(a,b), c(c, d)))$p.value #Calculate p-value
}

#Generate plot
dat.enrichment$pp<-factor(dat.enrichment$pp, levels=c(0.01, 0.05, 0.1,0.15,0.2, 0.25,0.3,0.4,0.5 ))
dat.enrichment<-subset(dat.enrichment, pp%in%c(0.01,0.05, 0.1,0.2))
ggplot(dat.enrichment, aes(fill=pp, y=fold, x=pp)) + 
  geom_bar(position="dodge", stat="identity", color="black")+
  theme_bw()+scale_fill_brewer(palette = "Blues")+ geom_hline(yintercept=1, linetype="dashed", color = "black")+
  xlab("PICS threshold")+ylab("Fold Enrichment")+theme(legend.position = "none")+
  geom_text(aes(label=a), position="dodge", vjust=1)+
  geom_text(aes(label=round(-log10(p),2)), position="dodge", vjust=-0.75)





