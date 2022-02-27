library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(GenomicRanges)
library(stringr)
library(data.table)

dir<-"/dir/T_cell_MPRA"


######################
####Run motifbreakR###
#Make input bed file
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra,project=="TGWAS" , select=c(chr, snp_start, snp_end, SNP, mpra_sig, LogSkew))
mpra$name<-paste(mpra$chr, mpra$snp_end, str_split_fixed(mpra$SNP, "\\:", 4)[,3], str_split_fixed(mpra$SNP, "\\:", 4)[,4],sep=":")
mpra$score<-0
mpra$strand<-"+"
mpra.bed<-subset(mpra, select=c(chr, snp_start, snp_end, name,score, strand))
write.table(mpra.bed, paste0(dir, "/tf/motifbreakr/mpra.motifbreakr.bed"), row.names = F, col.names = F, sep="\t", quote=F)

#import the BED file
snps.mb.frombed <- snps.from.file(file = paste0(dir, "/tf/motifbreakr/mpra.motifbreakr.bed"),
                                  dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                  format = "bed")

#Run motifbreakr and write results out
results <- motifbreakR(snpList =snps.mb.frombed, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-5,
                       method = "log")
dat.results<-data.frame(results, stringsAsFactors = F)
write.table(dat.results, paste0(dir, "/tf/motifbreakr/mpra.motifbreakr.method_log.p1e5.results.txt"), row.names = F, col.names = T, quote=F, sep="\t")


#####Evaluate overlap of MPRA variants with motifbreakR#######
#Read in main MPRA table
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")
mpra<-subset(mpra, project=="TGWAS", select=c( chr, snp_start, snp_end, SNP,A.log2FC,B.log2FC,LogSkew,mpra_sig))
mpra.se<-makeGRangesFromDataFrame(mpra, seqnames = "chr", start.field = "snp_start", end.field = "snp_end",keep.extra.columns=TRUE)

#Read in Motifbreakr data and merge with MPRA
dat.motifbreakr<-read.delim(paste0(dir, "/tf/motifbreakr/mpra.motifbreakr.method_log.p1e5.results.txt"),header=T, stringsAsFactors = F,sep="\t")
dat.motifbreakr<-subset(dat.motifbreakr,  select=c(seqnames, end, REF, ALT,SNP_id,  geneSymbol,scoreRef ,scoreAlt, alleleDiff))
dat.motifbreakr$SNP<-gsub("chr", "", dat.motifbreakr$SNP_id)
dat.motifbreakr<-merge(dat.motifbreakr, subset(mpra, select=c(SNP, A.log2FC,B.log2FC,LogSkew,mpra_sig)), by="SNP", all.x=T, all.y=F)

#Read in TF ChIP binding sites from ENCODE
dat.tf<-fread("https://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz")
names(dat.tf)[c(1:4,7)]<-c("chr", "start", "end", "tf", "cell")

#subset ENCODE TF data for txn factors identified in the motifbreakR calls
dat.tf<-subset(dat.tf, tf%in%dat.motifbreakr$geneSymbol, select=c(chr, start, end, tf, cell))
tf.se<-makeGRangesFromDataFrame(dat.tf, seqnames = "chr", start.field = "start", end.field = "end",keep.extra.columns=TRUE)

#Intersect MPRA variants with ENCODE TF data
dat.tf.overlap<-data.frame(subsetByOverlaps(tf.se, mpra.se, type="any"), stringsAsFactors = F)

#Find the cell types from ENCODE data where there is a TF binding site at each motifbreakR hit
dat.motifbreakr$overlap<-0
for(i in 1:nrow(dat.motifbreakr)){
  if(nrow(subset(dat.tf.overlap, tf==dat.motifbreakr[i,]$geneSymbol & seqnames==dat.motifbreakr[i,]$seqnames & start<=dat.motifbreakr[i,]$end & end>=dat.motifbreakr[i,]$end   ))){
    dat.motifbreakr[i,]$overlap<-subset(dat.tf.overlap, tf==dat.motifbreakr[i,]$geneSymbol & seqnames==dat.motifbreakr[i,]$seqnames & start<=dat.motifbreakr[i,]$end & end>=dat.motifbreakr[i,]$end   )$cell[1]
  }  
}

#Generate unique MPRA sites
dat.motifbreakr.unique<-subset(dat.motifbreakr,  overlap!=0)
dat.motifbreakr.unique<-dat.motifbreakr.unique[order(-abs(dat.motifbreakr.unique$alleleDiff)),]
dat.motifbreakr.unique<-dat.motifbreakr.unique[!duplicated(dat.motifbreakr.unique$SNP_id),]

#Subset MPRA file for SNPs evaluated by motifbreakR
snp.list<-snps.mb.frombed$SNP_id
mpra<-mpra[paste0("chr", mpra$SNP)%in%snps.mb.frombed$SNP_id,]

#Plot proportions of variants that have allelic TF binding
#Produces Figure S4B
dat.plot<-data.frame( mpra_sig=c("nEnhancer_nSkew", "Enhancer_nSkew", "Enhancer_Skew"), proportion=0,stringsAsFactors = F)
for(i in 1:nrow(dat.plot)){
  dat.plot[i,]$proportion<-nrow(subset(dat.motifbreakr.unique, mpra_sig==dat.plot[i,]$mpra_sig))/nrow(subset(mpra, mpra_sig==dat.plot[i,]$mpra_sig))
}
dat.plot$mpra_sig<-factor(dat.plot$mpra_sig, levels=c("nEnhancer_nSkew", "Enhancer_nSkew", "Enhancer_Skew"))

ggplot(dat.plot, aes(fill=mpra_sig, y=proportion, x=mpra_sig)) + 
  geom_bar(position="dodge", stat="identity", color="black")+
  theme_bw()+scale_fill_brewer(palette = "Blues")+coord_cartesian(ylim=c(0,0.05))+
  theme(axis.text.x = element_text( angle = 45, hjust=1))+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")

#Plot MPRA allelic skew against TF allelic difference from motifbreakR for emVars
dat.motifbreakr.emvar<-subset(dat.motifbreakr.unique, mpra_sig==c("Enhancer_Skew"))

plot(x=dat.motifbreakr.emvar$alleleDiff, y=dat.motifbreakr.emvar$LogSkew, pch=16, cex=1.2 ,col="forestgreen", xlim=c(-6,6), ylim=c(-1.5,1.5))
abline(v=0, lty=2, col="blue")
abline(h=0, lty=2, col="blue")
abline(lm(dat.motifbreakr.emvar$LogSkew~dat.motifbreakr.emvar$alleleDiff))

#Calculate statistical significance for correlation between MPRA alleliic skew against TF allelic difference
summary(lm(dat.motifbreakr.emvar$LogSkew~dat.motifbreakr.emvar$alleleDiff))

#Test significance of overlap between TFs with allelic activity and MPRA variants
mpra<-read.delim(paste0(dir, "/annotate_mpra/mpra_data_merge.txt"), header=T, stringsAsFactors = F, sep="\t")

#Perform for pCREs
prop.test(x = c(nrow(subset(dat.motifbreakr.unique , mpra_sig=="Enhancer_nSkew")), nrow(subset(dat.motifbreakr.unique, mpra_sig=="nEnhancer_nSkew"))), 
          n =c(nrow(subset(mpra, snp_end-snp_start==1 & mpra_sig=="Enhancer_nSkew")), nrow(subset(mpra, snp_end-snp_start==1 & mpra_sig=="nEnhancer_nSkew"))))

#Perform for emVars
prop.test(x = c(nrow(subset(dat.motifbreakr.unique , mpra_sig=="Enhancer_Skew")), nrow(subset(dat.motifbreakr.unique, mpra_sig=="nEnhancer_nSkew"))), 
          n =c(nrow(subset(mpra, snp_end-snp_start==1 & mpra_sig=="Enhancer_Skew")), nrow(subset(mpra, snp_end-snp_start==1 & mpra_sig=="nEnhancer_nSkew"))))
