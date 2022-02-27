library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(stringr)

sample<-"CD8_Na"
dir<-"/dir/T_cell_MPRA"

#Read in raw counts
dat<-read.delim(paste0(dir, "/rnaseq/brbseq/output.dge.umis.detailed.run3+4.DESeq.renamed.txt"), header = T, stringsAsFactors = F, sep="\t")
dat<-dat[!duplicated(dat$Gene_name),]
rownames(dat)<-dat$Gene_name
dat<-subset(dat, select=-c(Gene_id, Gene_name))
if(sample=="CD4_Na"){
  dat<-subset(dat, select=c(WT_forKO_CD4_Na_rep1, WT_forKO_CD4_Na_rep2, WT_forKO_CD4_Na_rep3, KO_CD4_Na_rep1, KO_CD4_Na_rep2, KO_CD4_Na_rep3))
}else if(sample=="CD8_Na"){
  dat<-subset(dat, select=c(WT_forKO_CD8_Na_rep1, WT_forKO_CD8_Na_rep2, WT_forKO_CD8_Na_rep3, KO_CD8_Na_rep1, KO_CD8_Na_rep2, KO_CD8_Na_rep3))
}
names(dat)<-gsub("_forKO", "", names(dat))


#Make annotations
dat.annot<-data.frame(sample=names(dat), condition=NA,  rep=NA, stringsAsFactors = F)
dat.annot$sample<-gsub("_forKO", "",dat.annot$sample)
dat.annot$condition<-str_split_fixed(dat.annot$sample, "_rep", 3)[,1]
dat.annot$rep<-gsub("rep", "", str_split_fixed(dat.annot$sample, "_rep", 3)[,2])
rownames(dat.annot)<-dat.annot[,1]
dat.annot<-subset(dat.annot, select=-c(sample))

#Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = dat.annot,
                              design = ~ condition)
vsd <- assay(vst(dds, blind=FALSE))
write.table(vsd, paste0(dir,"/rnaseq/brb_seq/",sample, ".deseq.vst.counts.txt"), row.names = T, col.names = T, sep="\t", quote=F)


#Make PCA plots
vsd<-vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_bw()+theme(legend.position = "bottom")


###Make normalized dotplots
dat.expr<-plotCounts(dds, gene="Bach2", intgroup="condition", returnData = T)
write.table(dat.expr, paste0(dir, "/rnaseq/brb_seq/",sample, ".normalized_counts.txt"), row.names = F, col.names = T, quote=F, sep="\t")
ggplot(dat.expr)+geom_dotplot(aes(x=condition, y=count, fill=condition),binaxis='y', stackdir='center')+theme_bw()+ theme(legend.position="none")

#run DESeq and make GSEA files
dds.results <- DESeq(dds)
diff <- results(dds.results, c("condition", paste0("KO_", sample), paste0("WT_", sample)), format = "DataFrame")
diff<-data.frame(diff, stringsAsFactors = F)
diff.out<-data.frame(gene=row.names(diff), stat=diff$stat)
diff.out<-diff.out[!is.na(diff.out$stat),]
write.table(diff.out, paste0(dir, "/rnaseq/gsea/KO_", sample, ".vs.WT_", sample, ".rnk"), row.names = F, col.names = F, sep="\t", quote=F)
write.table(diff, paste0(dir, "/rnaseq/gsea/KO_", sample, ".vs.WT_", sample, ".diff.txt"), row.names = T, col.names = T, sep="\t", quote=F)


###Plot heatmap for Bach2-18del mouse CD8 Naive BRB-seq, Figure 5E
#Read in normalized read counts
dat<-read.delim(paste0(dir, "/rnaseq/brb_seq/CD8_Na.deseq.vst.counts.txt"), row.names = 1, header=T, stringsAsFactors = F, sep="\t")
dat$gene<-row.names(dat)
dat<-subset(dat, select=c(gene, WT_CD8_Na_rep1 ,WT_CD8_Na_rep2 ,WT_CD8_Na_rep3,KO_CD8_Na_rep1,KO_CD8_Na_rep2,KO_CD8_Na_rep3))

#read in DE-seq results and subset for genes with differential expression pval < 0.05
dat.diff<-read.delim(paste0(dir, "/rnaseq/gsea/KO_CD8_Na.vs.WT_CD8_Na.diff.txt"), row.names = 1, header=T, stringsAsFactors = F, sep="\t")
dat.diff$gene<-row.names(dat.diff)
dat.diff<-subset(dat.diff, !is.na(stat) & padj<0.05,select=c(gene, stat))
dat<-merge(dat, dat.diff, by="gene", all.x=F, all.y=F)
dat<-dat[order(-dat$stat),]
row.names(dat)<-dat$gene
dat.annot<-subset(dat, select=c(stat))

#Row and column scale data
dat.scale<-data.frame(t(scale(t(subset(dat, select=-c(gene, stat))))), stringsAsFactors = F)
dat.scale[is.na(dat.scale)]<-0

#Plot heatmap
pheatmap(dat.scale, cluster_rows=T, cluster_cols=F ,annotation_row = dat.annot)
