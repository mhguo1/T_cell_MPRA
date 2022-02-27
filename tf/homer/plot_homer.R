library(ggplot2)
library(stringr)

dir<-"/dir/T_cell_MPRA"

#Plot Homer TF enrichment p-values
dat<-read.delim(paste0(dir, "/tf/homer/mpra_emvar_200/knownResults.txt"), header=T, stringsAsFactors = F, sep="\t")
dat$motif<-str_split_fixed(str_split_fixed(dat$Motif.Name, "\\/", 2)[,2], "\\-ChIP", 2)[,1]
dat<-dat[grepl("ChIP", dat$Motif.Name),]
dat$rank<-seq(1, nrow(dat))
plot(x=dat$rank, y=-dat$Log.P.value, cex=1.2, pch=16)

ggplot()+geom_point(data=dat,aes(x=rank,y=-Log.P.value),color="grey20", size=1.5, pch=16)+
   coord_cartesian(ylim=c(0,25), xlim=c(0,nrow(dat)))+ 
  theme_classic()+  geom_text(data=subset(dat, as.numeric(rank)<6),  aes(x=rank,y=-Log.P.value,label=motif),nudge_x = 10, nudge_y = 0.5, hjust=0)

