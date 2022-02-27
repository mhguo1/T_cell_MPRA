library(stringr)
library(ggplot2)

###Run PICS enrichment. Produces Figures 1C, 2B, S6A and S7
dir<-"/dir/T_cell_MPRA"
gwas.order<- c("Crohns","MS","Psoriasis", "RA","T1D","UC", "IBD" )
emvar_loci<-F #TRUE if doing only emvar loci. 

mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
names(mpra)<-gsub("_CS_", "_PP_", names(mpra))

mpra<-subset(mpra, project=="TGWAS", select=c(chr, snp_end, ld_snp, lead_snp, r2, rsid, Crohns_pval,Crohns_pics,Crohns_PP_running,MS_pval,MS_pics,MS_PP_running,
                            Psoriasis_pval,Psoriasis_pics,Psoriasis_PP_running,RA_pval,RA_pics,RA_PP_running,T1D_pval,T1D_pics,T1D_PP_running,
                            UC_pval,UC_pics,UC_PP_running,IBD_pval,IBD_pics,IBD_PP_running, dhs_Tcell_merged, dhs_all, mpra_sig))

#Remove bad SNPs where it doesn't reach 5E-8 association p-value in the GWAS. Also remove MHC region
bad_snps<-c("22:50966914:T:C","3:105558837:G:A", "12:9905851:A:C",
            "13:40745693:G:A","16:1073552:A:G","17:38775150:C:T",
            "17:44073889:A:G","18:12830538:G:A","2:100764087:T:G",
            "21:36488822:T:C","21:45621817:A:G","6:127457260:A:G",
            "6:130348257:C:T","7:116895163:G:A","7:51028987:T:A",
            "2:204592021:G:A", "14:75961511:C:T")
mpra<-subset(mpra,  !(chr=="chr6" & snp_end>29691116 & snp_end<33054976) & !(lead_snp%in%bad_snps))

#For each MPRA variant, find the disease with the strongest association and its associated PICS data
mpra$top_pval<-NA #Top GWAS p-value for the MPRA variant
mpra$top_disease<-NA #Disease corresponding to top GWAS p-value
mpra$top_PP_running<-NA #Cummulative sum of posterior probabilities for that variant
mpra$top_pics<-NA #PICS probability for that variant in the top GWAS
for(i in 1:nrow(mpra)){ #Run through each MPRA variant
  
  top_pval<-max(mpra[i,grepl("_pval",names(mpra))], na.rm=T) #Find the top GWAS p-value
  top_disease<-str_split_fixed(names(mpra)[which(mpra[i,]==top_pval)][1], "\\_", 2)[1] #Find the disease corresponding to the top GWAS p-value
  
  #Write out GWAS and PICS data for top GWAS p-value
  mpra[i,]$top_pval<-top_pval
  mpra[i,]$top_disease<-top_disease
  mpra[i,]$top_PP_running<-mpra[i,paste0(top_disease, "_PP_running")]
  mpra[i,]$top_pics<-mpra[i,paste0(top_disease, "_pics")]
}
mpra$top_pics<-as.numeric(mpra$top_pics)
mpra$top_PP_running<-as.numeric(mpra$top_PP_running)

dat.pics<-mpra

#find all GWAS loci with emVar SNPs
if(emvar_loci==T){
  emvar_lead_snps<-unique(subset(dat.pics, mpra_sig=="Enhancer_Skew")$lead_snp)
  dat.pics<-subset(dat.pics, lead_snp%in%emvar_lead_snps)
}

#Calculate MPRA/DHS enrichments in PICS fine-mapping
dat.enrichment<-data.frame(pics=rep(c(0.01, 0.05, 0.1,0.15,0.2, 0.25,0.3,0.4,0.5 ), times=3), 
                     disease=rep(rep(c("all"), times=9), times=3), 
                     mpra=rep(c("mpra", "dhs", "mpra_dhs"), each=9), 
                     a=0, b=0, c=0, d=0,fold=0, p=0, stringsAsFactors = F)

for(i in 1:nrow(dat.enrichment)){
  if(dat.enrichment[i,]$mpra=="mpra"){ #Calculate MPRA enrichments in PICS fine-mapping
    a<-nrow(subset(dat.pics, mpra_sig=="Enhancer_Skew" &  top_pics > dat.enrichment[i,]$pics)) #emVar SNP with PICS fine-mapped
    b<-nrow(subset(dat.pics, mpra_sig=="Enhancer_Skew" &  top_pics <= dat.enrichment[i,]$pics)) #emVar SNP, but PICS not fine-mapped
    c<-nrow(subset(dat.pics, mpra_sig!="Enhancer_Skew" &  top_pics > dat.enrichment[i,]$pics)) #Not emVar SNP, but PICS fine-mapped
    d<-nrow(subset(dat.pics, mpra_sig!="Enhancer_Skew" &  top_pics <= dat.enrichment[i,]$pics)) #Not emVar SNP, and not PICS fine-mapped
  }
  if(dat.enrichment[i,]$mpra=="dhs"){ #Calculate DHS enrichments in PICS fine-mapping
    a<-nrow(subset(dat.pics, dhs_Tcell_merged==1 &   top_pics > dat.enrichment[i,]$pics )) #DHS peak overlapping PICS fine-mapped SNP
    b<-nrow(subset(dat.pics, dhs_Tcell_merged==1  &  top_pics <= dat.enrichment[i,]$pics)) #DHS peak not overlapping PICS fine-mapped SNP
    c<-nrow(subset(dat.pics, dhs_Tcell_merged==0 &  top_pics > dat.enrichment[i,]$pics)) #Not overlapping DHS peak, but PICS fine-mapped SNP
    d<-nrow(subset(dat.pics, dhs_Tcell_merged==0 &  top_pics <= dat.enrichment[i,]$pics)) #Not overlapping DHS peak and not PICS fine-mapped
  }
  if(dat.enrichment[i,]$mpra=="mpra_dhs"){ #Calcualte MPRA+DHS enrichments in PICS fine-mapping
    a<-nrow(subset(dat.pics, mpra_sig=="Enhancer_Skew" & dhs_Tcell_merged==1 &   top_pics > dat.enrichment[i,]$pics )) #emVar, overlapping DHS peak and PICS fine-mapped
    b<-nrow(subset(dat.pics, mpra_sig=="Enhancer_Skew" & dhs_Tcell_merged==1  &  top_pics <= dat.enrichment[i,]$pics)) #emVar and overlapping DHS peak, but not PICS fine-mapped
    c<-nrow(subset(dat.pics, (mpra_sig!="Enhancer_Skew" | dhs_Tcell_merged==0) &  top_pics> dat.enrichment[i,]$pics)) #Either not emVar or not overlapping DHS peak, but PICS fine-mapped
    d<-nrow(subset(dat.pics, (mpra_sig!="Enhancer_Skew" | dhs_Tcell_merged==0)&  top_pics <= dat.enrichment[i,]$pics)) #Either not emVar or not overlapping DHS peak, and not PICS fine-mapped
  }
  
  #Write out data
  dat.enrichment[i,]$a<-a
  dat.enrichment[i,]$b<-b
  dat.enrichment[i,]$c<-c
  dat.enrichment[i,]$d<-d
  dat.enrichment[i,]$fold<-(a/(a+b))/(c/(c+d)) #Calculate fold enrichment
  dat.enrichment[i,]$p<-fisher.test(rbind(c(a,b), c(c, d)))$p.value #Calculate enrichment p-value
}

#Generate plot
dat.enrichment$pics<-factor(dat.enrichment$pics, levels=c(0.01, 0.05, 0.1,0.15,0.2, 0.25,0.3,0.4,0.5 ))
dat.enrichment<-subset(dat.enrichment, pics%in%c(0.01,0.05, 0.1,0.2, 0.3,0.4,0.5 ))
  ggplot(dat.enrichment, aes(fill=pics, y=fold, x=pics)) + 
    geom_bar(position="dodge", stat="identity", color="black")+facet_grid(~mpra)+
    theme_bw()+scale_fill_brewer(palette = "Blues")+ geom_hline(yintercept=1, linetype="dashed", color = "black")+
    xlab("PICS threshold")+ylab("Fold Enrichment")+theme(legend.position = "none")+
    geom_text(aes(label=a), position="dodge", vjust=1)+
    geom_text(aes(label=round(-log10(p),2)), position="dodge", vjust=-0.75)




###Sensitivity and specificity calculations###
dat.pics<-mpra
dhs_loci<-F #TRUE if calculation only for loci where a GWAS SNP overlaps a DHS peak
if(dhs_loci==T){  
  dat.pics<-subset(dat.pics, lead_snp%in%subset(dat.pics, dhs_Tcell_merged>0)$lead_snp)
}

#Calculate credible sets
dat.pics<-dat.pics[order(dat.pics$lead_snp, -dat.pics$top_pics),] 
dat.pics<-subset(dat.pics, select=c(ld_snp, lead_snp, r2, top_PP_running, top_pics, dhs_all, dhs_Tcell_merged, mpra_sig))
dat.pics$CS_80<-0 
dat.pics$CS_90<-0 
dat.pics$CS_95<-0 
for(i in 1:nrow(dat.pics)){
  top_pics<-max(subset(dat.pics, lead_snp==dat.pics[i,]$lead_snp)$top_pics)

  if(dat.pics[i,]$top_pics==top_pics){
    dat.pics[i,]$CS_80<-1
    dat.pics[i,]$CS_90<-1 
    dat.pics[i,]$CS_95<-1 
  }else{
    if(dat.pics[i,]$top_pics>=0.01){
      if(dat.pics[i,]$top_PP_running<=0.8){
        dat.pics[i,]$CS_80<-1
      }    
      if(dat.pics[i,]$top_PP_running<=0.9){
        dat.pics[i,]$CS_90<-1 
      }
      if(dat.pics[i,]$top_PP_running<=0.95){
        dat.pics[i,]$CS_95<-1 
      }
    }
  }
}



dat.sens_spec<-data.frame(cs=rep(c(80,90,95)), disease="all", mpra=rep(c("mpra"), each=3), 
                          a=0, b=0, c=0, d=0,sensitivity=0, specificity=0, stringsAsFactors = F)

for(i in 1:nrow(dat.sens_spec)){
  a<-0 #Instantiate variable for MPRA emVars that are in PICS credible set 
  b<-0 #variable for MPRA emVars that are not PICS credible set 
  c<-0 #variable for MPRA non-emVars that are in PICS credible set 
  d<-0 #variable for MPRA non-emVars that are not PICS credible set 
  a_list<-c()
  b_list<-c()
  c_list<-c()
  d_list<-c()
  for(l in unique(dat.pics$lead_snp)){ #Perform this calculation for each GWAS locus
    temp.dat.pics<-subset(dat.pics, lead_snp==l) #Subset SNPs for that locus
    total_list<-nrow(temp.dat.pics)
    cs_list<-nrow(temp.dat.pics[temp.dat.pics[,paste0("CS_", dat.sens_spec[i,]$cs)]==1,])
    non_cs_list<-nrow(temp.dat.pics[temp.dat.pics[,paste0("CS_", dat.sens_spec[i,]$cs)]==0,])
    max_pics<-max(temp.dat.pics$top_pics)
    if(max_pics>=0.05 & cs_list<=25 & total_list>=5 & total_list<100 & non_cs_list>0 ){
      dat.cs<-temp.dat.pics[temp.dat.pics[,paste0("CS_", dat.sens_spec[i,]$cs)]==1,]
      if(nrow(subset(dat.cs, mpra_sig=="Enhancer_Skew"))>0){
        a<-a+1
        a_list<-c(a_list, l)
      }else{
        c<-c+1
        c_list<-c(c_list, l)
      }
      
      dat.non_cs<-temp.dat.pics[temp.dat.pics[,paste0("CS_", dat.sens_spec[i,]$cs)]==0,]
      if(nrow(subset(dat.non_cs, mpra_sig=="Enhancer_Skew"))>0){
        b<-b+1
        b_list<-c(b_list, l)
      }else{
        d<-d+1
        d_list<-c(d_list, l)
      }
    }
  }
  dat.sens_spec[i,]$a<-a
  dat.sens_spec[i,]$b<-b
  dat.sens_spec[i,]$c<-c
  dat.sens_spec[i,]$d<-d
  dat.sens_spec[i,]$sensitivity<-a/(a+c) #Calculate sensitivity
  dat.sens_spec[i,]$specificity<-d/(b+d) #Calculate specificity
}
