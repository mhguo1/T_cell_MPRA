library(SummarizedExperiment)
library(ggplot2)
library(data.table )
library(stringr)
library(RColorBrewer)
library(ggrepel)

#Download DHS data from https://zenodo.org/record/3838751#.X69tfEJKg6U
#DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz
#DHS_Index_and_Vocabulary_metadata.tsv
#dat_bin_FDR01_hg19.txt.gz

dir<-"/dir/T_cell_MPRA"

#Read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra, project=="TGWAS", select=c(chr, snp_start, snp_end, mpra_sig))

#Make into a Granges file
mpra$strand<-"+"
mpra.se<-makeGRangesFromDataFrame(mpra,  seqnames = "chr", start.field = "snp_start", end.field = "snp_end",keep.extra.columns=TRUE)


#####################################
###Merge DHS peaks by cell type######
#####################################
#Read in DHS position file
dhs.pos<- read.delim(gzfile("DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz"), header=T, stringsAsFactors = F)
dhs.pos<-subset(dhs.pos, select=c(seqname, start, end))
dhs.pos$peak<-paste0(dhs.pos$seqname, ":", dhs.pos$start, "_", dhs.pos$end) #Make a peak name for easier indexing
dhs.pos.se<-makeGRangesFromDataFrame(dhs.pos, seqnames = "seqname", start.field = "start", end.field = "end",keep.extra.columns=TRUE) #Make into GRanges format

#Find DHS peaks that overlap with MPRA
overlap_peaks<-data.frame(subsetByOverlaps(dhs.pos.se, mpra.se, type="any"), stringsAsFactors = F)$peak #Find DHS peaks that overlap with MPRA SNPs
overlap_index<-which(dhs.pos$peak%in%overlap_peaks) #Extract row indices of overlapping peaks from DHS file
dhs.pos<-dhs.pos[overlap_index,] #Subset DHS peaks file by index of overlapping peaks
rm(mpra.se)
rm(dhs.pos.se)

#Read in DHS peak matrix
dhs.dat<-fread("dat_bin_FDR01_hg19.txt.gz")

#Subset DHS peak matrix for only peaks that overlap with MPRA
dhs.dat<-dhs.dat[overlap_index,] #Subset DHS matrix file by index of overlapping peaks
dhs.dat<-cbind(dhs.pos[,c(1:3)],dhs.dat)

#Read in DHS sample information
sample.dat<-read.delim("DHS_Index_and_Vocabulary_metadata.tsv", sep="\t", header=T, stringsAsFactors = F)
sample.dat<-subset(sample.dat, !is.na(library.order),select=c(library.order,Biosample.name, System,Subsystem,  Organ, Biosample.type, Biological.state))

#For each cell type, merge DHS peaks from that cell type. Add a column for each cell type, and list it's peak membership for any sample in that cell type
dhs.merge.dat<-dhs.dat[,c(1:3)]
for(c in unique(sample.dat$Biosample.name)){
  libraries<-paste0("V", subset(sample.dat, Biosample.name==c)$library.order)
  if(length(libraries)>1){
    dhs.merge.dat$newcol<-apply(dhs.dat[,c(libraries)],1, max)
  }else{
    dhs.merge.dat$newcol<-dhs.dat[,c(libraries)]
  }
  names(dhs.merge.dat)<-gsub("newcol", c, names(dhs.merge.dat))
}

#################################
###Run MPRA/pCRE Enrichments#####

#Read in sample information
sample.dat<-read.delim("DHS_Index_and_Vocabulary_metadata.tsv", sep="\t", header=T, stringsAsFactors = F)
sample.dat<-subset(sample.dat, !is.na(library.order),select=c(library.order,Biosample.name, System,Subsystem,  Organ, Biosample.type, Biological.state))

#Generate peak indices for MPRA SNPs for faster downstream processing
dhs.pos<-dhs.merge.dat[,c(1:3)]
dhs.pos$index<-seq(1:nrow(dhs.pos))
dhs.merge.dat$index<-seq(1:nrow(dhs.merge.dat))
mpra$peak_index<-NA
for(i in 1:nrow(mpra)){ #Add in DHS peak index for any peak overlapping MPRA SNP
  mpra[i,]$peak_index<-subset(dhs.pos, seqname==mpra[i,]$chr & start<=mpra[i,]$snp_end & end >=mpra[i,]$snp_end)$index[1]  
}
pcre_index<-subset(mpra, mpra_sig%in%c("Enhancer_nSkew", "Enhancer_Skew") & !is.na(peak_index))$peak_index #peak indices for MPRA pCRE SNPs only
non_pcre_index<-subset(mpra, mpra_sig=="nEnhancer_nSkew" & !is.na(peak_index))$peak_index #peak indices for MPRA non-pCRE SNPs

#Fix some sample names
sample.dat$Biosample.name<-gsub("786_O","X786_O", sample.dat$Biosample.name)
sample.dat$Biosample.name<-gsub("SK\\-N\\-DZ\\+RA\\s\\(3 days\\)","SK.N.DZ.RA..3.days.", sample.dat$Biosample.name)
sample.dat$Biosample.name<-gsub("Namalwa\\s\\+\\sSendai\\svirus_2h","Namalwa...Sendai.virus_2h", sample.dat$Biosample.name)
sample.dat$Biosample.name<-gsub( "MCF10a_ER_SRC_6h_Tam/Washout 18 h","MCF10a_ER_SRC_6h_Tam.Washout.18.h",sample.dat$Biosample.name)

#Calculate enrichments of MPRA pCREs in DHS peaks for each cell type
dat.plot<-data.frame(Biosample.name= unique(sample.dat$Biosample.name), System=NA, Subsystem=NA, pcre_fold=0, pcre_p=0, stringsAsFactors = F)
for(i in 1:nrow(dat.plot)){
  
  #Extract system for that cell type
  dat.plot[i,]$System<-subset(sample.dat, Biosample.name==dat.plot[i,]$Biosample.name)$System[1]
  
  #Extract subsystem for that cell type
  if(dat.plot[i,]$System!="Hematopoietic"){
    dat.plot[i,]$Subsystem<-dat.plot[i,]$System
  }else{
    dat.plot[i,]$Subsystem<-paste0("Heme_", subset(sample.dat, Biosample.name==dat.plot[i,]$Biosample.name)$Subsystem[1])
  }
  if(dat.plot[i,]$Subsystem=="Heme_"){dat.plot[i,]$Subsystem<-"Heme_Other"}
  if(dat.plot[i,]$Subsystem=="Fetal Life Support"){dat.plot[i,]$Subsystem<-"Fetal"}
  
  #Find the column in the dhs data that correspond to that cell type
  cell_col<-which(names(dhs.merge.dat)==dat.plot[i,]$Biosample.name)
  
  #find which peak_indexes for the peaks in that cell type
  peak_index<-dhs.merge.dat[dhs.merge.dat[,cell_col]==1,]$index

  #Calculate enrichment
  a<-length(intersect(peak_index, pcre_index)) #Number of pCRE in DHS peak
  b<-length(intersect(peak_index, non_pcre_index)) #Number of non-pCRE in DHS peak
  c<-nrow(subset(mpra, mpra_sig%in%c("Enhancer_nSkew", "Enhancer_Skew")))-a #Number of pCRE not in DHS
  d<-nrow(subset(mpra, mpra_sig=="nEnhancer_nSkew"))-b #Number of non-pCRE not in DHS
  dat.plot[i,]$pcre_fold<-(a/(a+c))/(b/(b+d)) #Calcualte fold enrichment
  dat.plot[i,]$pcre_p<-fisher.test(rbind(c(a,b), c(c, d)), alternative="two.sided")$p.value #Calculate p-value
}

#Rank data by p-value
dat.plot<-dat.plot[order(dat.plot$pcre_p),]
dat.plot$rank<-seq(1:nrow(dat.plot))

#Generate plot
subsystem.order<-c("Cardiovascular","Connective","Digestive","Embryonic","Endocrine","Epithelial","Fetal","Genitourinary", "Hepatic","Integumentary","Musculoskeletal","Nervous","Renal","Respiratory","Heme_Erythroid","Heme_Myeloid","Heme_Other","Heme_Progenitor","Heme_Lymphoid","Heme_T-cell")
dat.plot$Subsystem<-factor(dat.plot$Subsystem, levels=subsystem.order)
dat.plot<-dat.plot[order(dat.plot$pcre_p),]
dat.plot$rank<-seq(1:nrow(dat.plot))
ggplot()+geom_point(data=dat.plot[!(grepl("Heme", dat.plot$Subsystem)),],aes(x=rank,y=-log10(pcre_p)),color="grey", size=1.5, pch=16)+
  geom_point(data=dat.plot[(grepl("Heme", dat.plot$Subsystem)),],aes(x=rank,y=-log10(pcre_p), color=Subsystem), size=3, pch=16)+
  coord_cartesian(ylim=c(0,40), xlim=c(0,nrow(dat.plot)))+ 
  scale_color_manual(values=brewer.pal(6, "Dark2"))+
  theme_classic()+  geom_text(data=subset(dat.plot, as.numeric(rank)<6),  aes(x=rank,y=-log10(pcre_p),label=Biosample.name),nudge_x = 10, nudge_y = 0.5, hjust=0)


