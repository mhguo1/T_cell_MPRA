library(pheatmap)



dir<-"/dir/rnaseq/"
###Draw heatmap of Tscm data based on diff in our data
dat<-read.delim(paste0(dir, "/GSE152379/RNA_Bach2gRNA.vs.pMK9_P14_Cl13_D7_StemLike.deseq.vst.counts.txt"), header=T, stringsAsFactors = F)
dat$gene<-row.names(dat)

#read in differentially expressed genes and merge
dat.diff<-read.delim(paste0(dir, "/gsea/KO_CD8_Na.vs.WT_CD8_Na.diff.txt"), row.names = 1, header=T, stringsAsFactors = F, sep="\t")
dat.diff$gene<-row.names(dat.diff)
dat.diff<-subset(dat.diff, !is.na(stat) & padj<0.05,select=c(gene, stat))
dat<-merge(dat, dat.diff, by="gene", all.x=F, all.y=F)
dat<-dat[order(-dat$stat),]
row.names(dat)<-dat$gene
dat.annot<-subset(dat, select=c(stat))


#Expression of genes in Tscms from empty vector vs Bach2 sgRNA for differentially upregulated genes in Bach218del vs. WT naïve CD8 T cells, 
dat.up<-subset(dat, stat>0)
dat.annot.up<-subset(dat.up, select=c( stat))
row.names(dat.annot.up)<-dat.up$gene
dat.annot.up[nrow(dat.annot.up)+1,1]<-mean(dat.annot.up$stat)
dat.scale.up<-data.frame(t(scale(t(subset(dat.up, select=-c(gene, stat))))), stringsAsFactors = F)
dat.scale.up[nrow(dat.scale.up)+1,]<-c( -2,2,1,1,-2,2)
dat.scale.up[is.na(dat.scale.up)]<-0
pheatmap(dat.scale.up, cluster_rows=T, cluster_cols=F,annotation_row = dat.annot.up, cutree_rows = 2)


#Expression of genes in Tscms from empty vector vs Bach2 sgRNA for differentially downregulated genes in Bach218del vs. WT naïve CD8 T cells
dat.down<-subset(dat, stat<0)
dat.annot.down<-subset(dat.down, select=c( stat))
row.names(dat.annot.down)<-dat.down$gene
dat.annot.down[nrow(dat.annot.down)+1,1]<-mean(dat.annot.down$stat)
dat.scale.down<-data.frame(t(scale(t(subset(dat.down, select=-c(gene, stat))))), stringsAsFactors = F)
dat.scale.down[nrow(dat.scale.down)+1,]<-c( -2,2,1,1,-2,2)
dat.scale.down[is.na(dat.scale.down)]<-0
ann_colors = list(stat = c("white", "firebrick"))
pheatmap(dat.scale.down, cluster_rows=T, cluster_cols=F,annotation_row = dat.annot.down, annotation_colors=ann_colors, cutree_rows = 2)
