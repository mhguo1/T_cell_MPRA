library(stringr)
library(IRanges)
library(SummarizedExperiment)
library(ggplot2)
library(pheatmap)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

dir<-"/dir/T_cell_MPRA"

#########################################################
####Make Histogram of distance of MPRA variant to TSS###

#read in MPRA data
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra, project=="TGWAS" & !is.na(LogSkew))

#Expand out MPRA data
mpra.all<-mpra
mpra.all$mpra_sig_expanded<-"all"
mpra.pcre<-subset(mpra, mpra_sig%in%c("Enhancer_nSkew","Enhancer_Skew"))
mpra.pcre$mpra_sig_expanded<-"pcre"
mpra.emvar<-subset(mpra, mpra_sig=="Enhancer_Skew")
mpra.emvar$mpra_sig_expanded<-"emvar"
mpra.expanded<-rbind(mpra.all, mpra.pcre, mpra.emvar)

#Annotate MPRA functional annotations
mpra.all.se<-makeGRangesFromDataFrame(mpra.expanded,  seqnames = "chr", start.field = "snp_start", end.field = "snp_end",keep.extra.columns=TRUE)
peakAnno <- as.data.frame(annotatePeak(mpra.all.se, tssRegion=c(-3000, 3000),
                                       TxDb=txdb, annoDb="org.Hs.eg.db"))



#Plot distance to TSS as histogram
dat.plot<-data.frame(start=c(0, 1000, 3000, 10000, 20000 ),
                     end=c(1000, 3000, 10000,  20000,1000000), 
                     all=0, pcre=0, emvar=0, stringsAsFactors = F)
for(i in 1:nrow(dat.plot)){
  dat.plot[i,]$all<-nrow(subset(peakAnno, distanceToTSS>dat.plot[i,]$start & distanceToTSS<=dat.plot[i,]$end ))/nrow(peakAnno)
  dat.plot[i,]$pcre<-nrow(subset(peakAnno, mpra_sig%in%c("Enhancer_nSkew", "Enhancer_Skew") & distanceToTSS>dat.plot[i,]$start & distanceToTSS<=dat.plot[i,]$end ))/nrow(subset(peakAnno, mpra_sig%in%c("Enhancer_nSkew", "Enhancer_Skew")))
  dat.plot[i,]$emvar<-nrow(subset(peakAnno, mpra_sig%in%c("Enhancer_Skew") & distanceToTSS>dat.plot[i,]$start & distanceToTSS<=dat.plot[i,]$end ))/nrow(subset(peakAnno, mpra_sig%in%c("Enhancer_Skew")))
}

dat.plot$label=paste(dat.plot$start, dat.plot$end, sep="_")
dat.plot<-subset(dat.plot, select=-c(start, end))
dat.plot<-melt(dat.plot, id.vars="label")
dat.plot$label<-factor(dat.plot$label, levels=dat.plot$label[1:5])


ggplot(dat.plot, aes(x=label, y=value, fill=variable))+ 
  geom_bar(stat="identity", width=0.9, position = "dodge") +
  theme_bw()+scale_fill_brewer(palette = "Blues")+
  theme(axis.text.x = element_text( angle = 45, hjust=1))+
  xlab("Distance to TSS (bp)") + ylab("Proportion of variants")


