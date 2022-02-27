library(readxl)
library(ggplot2)
library(data.table)
library(stringr)
library(IRanges)
library(SummarizedExperiment)

dir<-"/dir/T_cell_MPRA"

#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")

mpra<-subset(mpra,project=="TGWAS", select=c(chr, snp_start, snp_end, ID, mpra_sig, dhs_Tcell_merged))
mpra.se<-makeGRangesFromDataFrame(mpra, seqnames = "chr", start.field = "snp_start", end.field = "snp_end",keep.extra.columns=TRUE)

#####################################
###Add in histone mark annotations###

#Read in manifest file containing file paths to ENCODE histone ChIP bed files for T cells
manifest.histone<-data.frame(read_excel("histone_manifest.xlsx", sheet=1), stringsAsFactors = F)

#Add histone data to MPRA dataframe
histone.list<-unique(manifest.histone$Histone)
for(h in histone.list){ #Run for each histone mark
  bed_url.list<-subset(manifest.histone, Histone==h )$link #Grab urls for all bed files for that histone mark
  
  bed<-data.frame(V1=character(), V2=numeric(), V3=numeric(), stringsAsFactors = F) #Create empty bed file
  for(url in bed_url.list){ #Read in and concatenate each bed file of the histone mark
    temp<-data.frame(fread(url), stringsAsFactors = F)
    temp<-temp[temp$V1%in%paste0("chr", seq(1,22)),c(1:3)]
    bed<-rbind(temp, bed)
  }
  names(bed)<-c("chr", "start", "end")
  bed.se<-makeGRangesFromDataFrame(bed, seqnames = "chr", start.field = "start", end.field = "end") #Create GRange for concatenated bed file
  mpra$newcol<-ifelse(mpra$ID%in%(subsetByOverlaps( mpra.se, bed.se,type="any")$ID),1,0) #Insert new column in MPRA file with whether MPRA variant overlaps the histone mark
  names(mpra)[ncol(mpra)]<-h
}

#####################################
###Add in CAGE annotations###
cage<-data.frame(fread("https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz")) #read in CAGE data from FANTOM consortium
cage.se<-makeGRangesFromDataFrame(cage, seqnames = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = F) #Make GRange for CAGE data
mpra$cage<-ifelse(mpra$ID%in%(subsetByOverlaps( mpra.se, cage.se,type="any")$ID),1,0) #Add in new column in MPRA file with whether or not it is a CAGE-based enhancer

#Plot proportion of MPRA variants overlapping chromHMM, DHS, or CAGE
mpra$dhs<-mpra$dhs_Tcell_merged #Get merged T cell DHS column

dat.enrichment<-data.frame(mark=rep(c(histone.list, "cage", "dhs"), each=2), mode=rep(c("pcre", "emvar"), times=8), a=0, b=0, c=0, d=0, fold=0, p=0,stringsAsFactors = F)

for(i in 1:nrow(dat.enrichment)){
  cell_col<-which(names(mpra)==dat.enrichment[i,]$mark) #Get column of MPRA file corresponding to current anntotation
  if(dat.enrichment[i,]$mode=="emvar"){ #Run enrichment for emVars in the annotation
    a<-nrow(mpra[mpra$mpra_sig=="Enhancer_Skew" & mpra[,cell_col]==1,])
    b<-nrow(mpra[mpra$mpra_sig=="Enhancer_Skew"& mpra[,cell_col]==0,])
    c<-nrow(mpra[mpra$mpra_sig!="Enhancer_Skew" & mpra[,cell_col]==1,])
    d<-nrow(mpra[mpra$mpra_sig!="Enhancer_Skew" & mpra[,cell_col]==0,])
  }
  if(dat.enrichment[i,]$mode=="pcre"){ #Run enrichment for pCREs in the annotation
    a<-nrow(mpra[mpra$mpra_sig%in%c("Enhancer_Skew", "Enhancer_nSkew") & mpra[,cell_col]==1,])
    b<-nrow(mpra[mpra$mpra_sig%in%c("Enhancer_Skew", "Enhancer_nSkew")& mpra[,cell_col]==0,])
    c<-nrow(mpra[mpra$mpra_sig=="nEnhancer_nSkew" & mpra[,cell_col]==1,])
    d<-nrow(mpra[mpra$mpra_sig=="nEnhancer_nSkew" & mpra[,cell_col]==0,])
  }
  dat.enrichment[i,]$a<-a
  dat.enrichment[i,]$b<-b
  dat.enrichment[i,]$c<-c
  dat.enrichment[i,]$d<-d
  dat.enrichment[i,]$fold<-(a/(a+c))/(b/(b+d)) #Calculate enrichment
  dat.enrichment[i,]$p<-fisher.test(rbind(c(a,b), c(c, d)), alternative="two.sided")$p.value #Calculate p-value
}
write.table(dat.enrichment, "ENCODE_histone_enrichment.txt", row.names = F, col.names = T, sep="\t", quote=F)
dat.enrichment$mode<-factor(dat.enrichment$mode, levels=c("pcre", "emvar"))
dat.enrichment$mark<-factor(dat.enrichment$mark, levels=c("H3K4me1","H3K4me3", "H3K9me3", "H3K27ac","H3K27me3", "H3K36me3", "cage", "dhs"))
dat.enrichment$significant<-ifelse(dat.enrichment$p<(0.05/8), "*", "")

#Paired barplot of histone overlaps
ggplot(dat.enrichment, aes(fill=mode, y=fold, x=mark)) + 
  geom_bar(position="dodge", stat="identity", color="black")+
  theme_bw()+scale_fill_brewer(palette = "Paired")+coord_cartesian(ylim=c(0,4))+
  theme(axis.text.x = element_text( angle = 45, hjust=1))+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  geom_text(aes(label=significant), position=position_dodge(width=0.9), vjust=-0.25)



##############################
###Overlap with chromHMM#####

#Read in chromHMM file from ENCODE
chromhmm<-fread("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/E034_18_core_K27ac_dense.bed.gz", data.table = F, stringsAsFactors = F, header=F)
chromhmm<-chromhmm[,c(1:4)]
names(chromhmm)<-c("chr", "start", "end", "chromhmm")

#Add in membership of chromHMM segmentations for MPRA variants
mpra$chromhmm<-NA
for(i in 1:nrow(mpra)){
  mpra[i,]$chromhmm<-subset(chromhmm, chr==mpra[i,]$chr & start<mpra[i,]$snp_end & end>=mpra[i,]$snp_end)$chromhmm[1]
}

#Calculate enrichments of pCREs and emVars in chromHMM segmentations
dat.enrichment<-data.frame(chromhmm=rep(unique(chromhmm$chromhmm), 2),mode=rep(c("pcre", "emvar"), each=length(unique(chromhmm$chromhmm))), a=0, b=0, c=0, d=0,fold=0, p=0, stringsAsFactors = F)
for(i in 1:nrow(dat.enrichment)){
  if(dat.enrichment[i,]$mode=="pcre"){ #Run enrichment for pCREs in the chromHMM segmentation
    a<-nrow(subset(mpra, mpra_sig!="nEnhancer_nSkew" & chromhmm==dat.enrichment[i,]$chromhmm))
    b<-nrow(subset(mpra, mpra_sig!="nEnhancer_nSkew" & chromhmm=="18_Quies"))
    c<-nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew" & chromhmm==dat.enrichment[i,]$chromhmm))
    d<-nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew" & chromhmm=="18_Quies"))
  }
  if(dat.enrichment[i,]$mode=="emvar"){ #Run enrichment for emVars in the chromHMM segmentation
    a<-nrow(subset(mpra, mpra_sig=="Enhancer_Skew" & chromhmm==dat.enrichment[i,]$chromhmm))
    b<-nrow(subset(mpra, mpra_sig=="Enhancer_Skew" & chromhmm=="18_Quies"))
    c<-nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew" & chromhmm==dat.enrichment[i,]$chromhmm))
    d<-nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew" & chromhmm=="18_Quies"))
  }
  dat.enrichment[i,]$a<-a
  dat.enrichment[i,]$b<-b
  dat.enrichment[i,]$c<-c
  dat.enrichment[i,]$d<-d
  dat.enrichment[i,]$fold<-(a/(a+b))/(c/(c+d)) #Calculate enrichment
  dat.enrichment[i,]$p<-fisher.test(rbind(c(a,b), c(c, d)), alternative="two.sided")$p.value #calculate p-value
}
write.table(dat.enrichment, "ENCODE_chromhmm_enrichment.txt", row.names = F, col.names = T, sep="\t", quote=F)

#Generate plots
dat.enrichment$segmentation<-as.numeric(str_split_fixed(dat.enrichment$chromhmm, "\\_", 2)[,1])
dat.enrichment<-dat.enrichment[order(dat.enrichment$segmentation),]
dat.enrichment$chromhmm<-factor(dat.enrichment$chromhmm, levels=unique(dat.enrichment$chromhmm))
dat.enrichment$mode<-factor(dat.enrichment$mode, levels=c("pcre", "emvar"))
dat.enrichment$significant<-ifelse(dat.enrichment$p<(0.05/36), "***", "")
dat.enrichment[dat.enrichment$p<0.05 & dat.enrichment$p>=0.05/36,]$significant<-"*"

#Paired barplot of chromHMM
ggplot(dat.enrichment, aes(fill=mode, y=fold, x=chromhmm)) + 
  geom_bar(position="dodge", stat="identity", color="black")+
  theme_bw()+scale_fill_brewer(palette = "Paired")+coord_cartesian(ylim=c(0,20))+
  theme(axis.text.x = element_text( angle = 45, hjust=1))+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  geom_text(aes(label=significant), position=position_dodge(width=0.9), vjust=-0.25)
