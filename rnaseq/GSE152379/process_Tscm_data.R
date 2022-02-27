library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(stringr)
require(biomaRt)
require(readxl)
dir<-"/dir/T_cell_MPRA"

convertMouseGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  print(head(humanx))
  return(humanx)
}

#Read in RNA-seq count matrix for GSE152379
dat.counts<-read.delim(paste0(dir, "/rnaseq/GSE152379/GSE152379.kallisto.merged.txt"), header=T, stringsAsFactors = F, sep="\t")
dat.genes<-subset(dat.counts, select=c(target_id,ensg, gene))
dat.counts<-subset(dat.counts, select=c(target_id,ensg, gene,RNA_pMK9_P14_Cl13_D7_StemLike_1, RNA_pMK9_P14_Cl13_D7_StemLike_2, RNA_pMK9_P14_Cl13_D7_StemLike_3,RNA_Bach2gRNA_P14_Cl13_D7_StemLike_1, RNA_Bach2gRNA_P14_Cl13_D7_StemLike_2, RNA_Bach2gRNA_P14_Cl13_D7_StemLike_3))

#Aggregate counts across transcripts of the same gene
dat<-data.frame(gene=unique(dat.counts$gene), RNA_pMK9_P14_Cl13_D7_StemLike_1=0, RNA_pMK9_P14_Cl13_D7_StemLike_2=0, RNA_pMK9_P14_Cl13_D7_StemLike_3=0,RNA_Bach2gRNA_P14_Cl13_D7_StemLike_1=0, RNA_Bach2gRNA_P14_Cl13_D7_StemLike_2=0, RNA_Bach2gRNA_P14_Cl13_D7_StemLike_3=0,stringsAsFactors = F)
for(i in 1:nrow(dat)){
  dat[i,c(2:7)]<-colSums(dat.counts[dat.counts$gene==dat[i,]$gene, c(4:9)])
}
rownames(dat)<-dat$gene
dat<-subset(dat, select=-c(gene))
dat[is.na(dat)]<-0
dat<-round(dat)


#Make annotation file for DESeq
dat.annot<-data.frame(sample=names(dat), condition=c(rep("RNA_pMK9_P14_Cl13_D7_StemLike",3),rep("RNA_Bach2gRNA_P14_Cl13_D7_StemLike",3)) , rep=c(1,2,3,1,2,3), stringsAsFactors = F)
rownames(dat.annot)<-dat.annot[,1]
dat.annot<-subset(dat.annot, select=-c(sample))

#Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = dat.annot,
                              design = ~ condition)
vsd <- assay(vst(dds, blind=FALSE)) #vst normalize RNA-seq data
write.table(vsd, paste0(dir, "/rnaseq/GSE152379/RNA_Bach2gRNA.vs.pMK9_P14_Cl13_D7_StemLike.deseq.vst.counts.txt"), row.names = T, col.names = T, sep="\t", quote=F)


#run DESeq
dds.results <- DESeq(dds)
diff <- results(dds.results, c("condition", "RNA_Bach2gRNA_P14_Cl13_D7_StemLike", "RNA_pMK9_P14_Cl13_D7_StemLike"), format = "DataFrame")
diff<-data.frame(diff, stringsAsFactors = F)
write.table(diff, paste0(dir, "/rnaseq/GSE152379/RNA_Bach2gRNA.vs.pMK9_P14_Cl13_D7_StemLike.diff.txt"), row.names = T, col.names = T, sep="\t", quote=F)

#Write out list of top 200 significantly differentially upregulated genes upon Bach2 knockdown
diff.up<-(subset(diff, stat>0 & padj<0.05))
diff.up<-diff.up[order(-diff.up$stat),]
musgenes<-row.names(diff.up[1:200,])
humgenes<-convertMouseGeneList(musgenes)
write.table(humgenes, paste0(dir, "/rnaseq/GSE152379/RNA_Bach2gRNA.vs.pMK9_P14_Cl13_D7_StemLike.diff_genes.txt"), row.names = F, col.names = F, sep="\t", quote=F)

#Write out list of top 200 significantly differentially downregulated genes upon Bach2 knockdown
diff.down<-(subset(diff, stat<0 & padj<0.05))
diff.down<-diff.down[order(diff.down$stat),]
musgenes<-row.names(diff.down[1:200,])

humgenes<-convertMouseGeneList(musgenes)
write.table(humgenes, paste0(dir, "/rnaseq/GSE152379/pMK9_P14_Cl13_D7_StemLike.vs.RNA_Bach2gRNA.diff_genes.txt"), row.names = F, col.names = F, sep="\t", quote=F)
