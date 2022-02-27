library(stringr)
library(data.table)

#This script takes the GWAS catalog and extracts SNPs that were tested in the MPRA.

linker<-fread("all.rsid.linker",header=F, stringsAsFactors = F, data.table = F) 
#This linker file can be make by concatenating all the 1000 Genomes vcf files and extracting their rsIDs and SNPIDs (in CHROM:POS:REF:ALT format)
#zcat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v ^# | awk '{print $3"\t"$1":"$2":"$4":"$5}' | sort | uniq > all.rsid.linker


names(linker)<-c("SNPS", "SNPID")
linker$chr_pos<-paste0(str_split_fixed(linker$SNPID, "\\:", 4)[,1], ":",str_split_fixed(linker$SNPID, "\\:", 4)[,2]) #Create a temporary field in linker file with just chromosome and position

#Run through each trait and select out GWAS catalog SNPs that were tested in the MPRA
for(d in c("MS", "Crohns", "IBD", "UC", "Celiac", "SLE", "Psoriasis", "PBC", "RA", "T1D")){
  if(d=="MS"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0003929-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Multiple sclerosis" | DISEASE.TRAIT=="Marginal zone lymphoma or multiple sclerosis")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
    
    #Some SNPs were tested in MPRA and associated with MS but are missing in 1000 Genomes. Use a nearby SNP in perfect LD to circumvent this issue 
    dat$SNPS<-gsub("rs12832171", "rs767455", dat$SNPS) 
    dat$SNPS<-gsub("chr11:14868316", "rs116970203", dat$SNPS) 
    dat$SNPS<-gsub("rs4820955","rs2027982", dat$SNPS)
  }
  if(d=="Crohns"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0000384-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Crohn's disease")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  }
  if(d=="IBD"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0000384-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Inflammatory bowel disease" | DISEASE.TRAIT=="Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
    dat<-rbind(dat, c("rs55722650", 4.61E-41)) #Manually add this SNP to GWAS catalog, from PMID 26192919
  }
  if(d=="UC"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0000384-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Ulcerative colitis")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  }
  if(d=="Celiac"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0001060.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Celiac disease")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  }
  if(d=="SLE"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0004537-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Systemic lupus erythematosus")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  }
  if(d=="T1D"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_0009706-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Type 1 diabetes")
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  }
  if(d=="PBC"){
    dat<-read.delim("gwas-association-downloaded_2020-08-10-EFO_1001486.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Primary biliary cirrhosis" | DISEASE.TRAIT=="Primary biliary cholangitis" )
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  }
  if(d=="RA"){
    dat<-read.delim("gwas-association-downloaded_2020-10-07-EFO_0009460-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Rheumatoid arthritis" )
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
  } 
  if(d=="Psoriasis"){
    dat<-read.delim("gwas-association-downloaded_2020-11-08-EFO_1001494-withChildTraits.tsv", header=T, stringsAsFactors = F)
    dat<-subset(dat, DISEASE.TRAIT=="Psoriasis" )
    dat<-subset(dat, select=c("SNPS", "P.VALUE"))
    dat<-dat[!grepl("x", dat$SNPS),]
  } 
  
  #Fix issues with some SNPs in catalog not having an rsID, but rather being "chr:pos"
  for(i in 1:nrow(dat)){
    if(grepl("chr", dat[i,]$SNPS)){
      dat[i,]$SNPS<-gsub("\\ ", "",dat[i,]$SNPS )
      dat[i,]$SNPS<-gsub("chr", "",dat[i,]$SNPS )
      if(dat[i,]$SNPS%in%linker$chr_pos){
        dat[i,]$SNPS<-subset(linker, chr_pos==dat[i,]$SNPS)$SNPS[1] 
      }
    }
  }
  
  #Manually fix a few catalog SNPs where the rsID has changed in 1000 Genomes
  dat$SNPS<-gsub("rs114934997","rs112768831", dat$SNPS ) 
  dat$SNPS<-gsub("rs917911", "ss1388113193" , dat$SNPS )  
  dat$SNPS<-gsub("rs11052877", "ss1388113192" , dat$SNPS ) 
  dat$SNPS<-gsub("rs7977720", "ss1388111988" , dat$SNPS )  
  dat$SNPS<-gsub("rs76901167", "rs9271100", dat$SNPS) 
  dat$SNPS<-gsub("rs67934705", "rs112741635", dat$SNPS) 
  dat$SNPS<-gsub("rs34723276", "rs74449127", dat$SNPS)
  
  #replace catalog SNP with a nearby SNP in perfect LD
  dat$SNPS<-gsub("rs12663356", "rs12661871", dat$SNPS) 
  
  #Calculate -log10(p) and cap at 350
  dat$P.VALUE<-(-log10(as.numeric(dat$P.VALUE)))
  dat$P.VALUE[is.infinite(dat$P.VALUE)]<-350

  #Remove duplicate SNPS
  dat<-dat[order(-dat$P.VALUE),]
  dat<-dat[!(duplicated(dat$SNPS)),]
  
  #Merge with linker file to get 1KG SNPID from the catalog rsIDs
  dat<-merge(dat, linker, by="SNPS", all.x=F, all.y=F)
  dat$SNPID<-paste0("chr", dat$SNPID)
  dat<-dat[!(grepl("chrX", dat$SNPID)),]
  
  #Output data
  dat<-dat[,c(3,2,1)]
  write.table(dat, paste0( d, ".catalog_sumstats.txt"), row.names = F, col.names = F, sep=" ", quote=F)
}
