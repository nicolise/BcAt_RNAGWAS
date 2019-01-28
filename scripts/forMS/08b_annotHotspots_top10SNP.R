#08_annotHotspots.R
#Nicole E Soltis
#01/24/19
#-------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS")
TopSNPAll <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/Top10SNP_BcAt_20topSNPs_HotSpot_rmPermutSNP.csv")
plotSNP <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/Top10SNP_BcAt_allhotspots_HotSpot_rmPermutSNP.csv")
# make stringent TopSNP with At > 150, Bc > 20

#split out chr_snp
library(stringr)
blah <- as.data.frame(str_split_fixed(TopSNPAll$chr_snp, "_", 2))
TopSNPAll <- cbind(TopSNPAll, blah)
names(TopSNPAll)[7] <- "chr"
names(TopSNPAll)[8] <- "pos"
TopSNPAll$chr <- as.numeric(paste(TopSNPAll$chr))
TopSNPAll$pos <- as.numeric(paste(TopSNPAll$pos))
TopSNPAll <- TopSNPAll[order(c(TopSNPAll$chr, TopSNPAll$pos)),]
TopSNPAll <- TopSNPAll[1:45,]

#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]

#calculate gene center
#calculate distance gene center to SNP
#add gene with min distance
#range +-1 kb around each snp: lowrange toprange
#match snp chromosome.id to gene V1

#associate each plant SNP with nearest gene from my.gtf (this is B05.10 gene annotation)
for (j in c(1:6,8:14,16)){ #none on chr 7,15,17,18
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  TopSNP.sub <- TopSNPAll[TopSNPAll$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (i in c(1:nrow(TopSNP.sub))){
    this.snp <- as.numeric(TopSNP.sub[i,8]) #check that this is the positional info/ SNP location
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(TopSNP.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, TopSNP.genes <- all.genes, TopSNP.genes <- rbind(TopSNP.genes, all.genes))
}
#now only keep genes if nearest end is within 1kb of SNP
TopSNP.genes$pos <- as.numeric(unlist(TopSNP.genes$pos))
TopSNP.genes$closest.end <- pmin(abs(TopSNP.genes$pos - TopSNP.genes$V4),abs(TopSNP.genes$pos - TopSNP.genes$V5)) 
TopSNP.genes.sub <- TopSNP.genes[TopSNP.genes$closest.end < 2000,]
#at this point, includes multiple transcripts per gene and multiple SNPs per gene
setwd("~/Projects/BcAt_RNAGWAS")
write.csv(TopSNP.genes, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/BcAt_top10SNP_hotspot_Genes.csv")
blah <- as.data.frame(unique(TopSNP.genes$V12))
write.csv(blah, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/BcAt_top10SNP_hotspot_geneNames.csv")
#--------------------------------------------------------------------------------
#look for function annotations for these genes!
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
full.snp.genes <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/BcAt_top10SNP_hotspot_Genes_ed.csv")
funclist <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/BcAt_top10SNP_hotspot_geneNames_botport.csv")
funclist <- funclist[,-c(1)]
names(funclist)[1] <- "Gene"
#could filter out ind-SNP info at this stage and only keep one record per gene... also sum across SNPs for hotspot info
#and take minimum of closest.end
full.genes <- full.snp.genes
full.genes <- unique(full.genes)
names(full.genes)
library(dplyr)
#50 unique genes
fg.sum <- full.genes %>%
  group_by(Gene) %>% 
  summarise(Gene.B = sum(Gene.B), 
            Gene.A = sum(Gene.A),
            V4 = min(V4),
            v5 = max(V5),
            midgene = mean(midgene),
            closest.end = min(closest.end)
            )
fg.deets <- full.genes[,c("Gene","chr","Cat","V14")]
fg.deets <- unique(fg.deets)
fg.sum <- merge(fg.sum, fg.deets, by="Gene")
fg.sum <- unique(fg.sum)

full.gene.funcs <- merge(fg.sum, funclist, by="Gene", all.x=T)
write.csv(full.gene.funcs,"data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/BcAt_top10SNP_hotspots_funcannot.csv")

##START HERE
#and try connecting this list to Vivian's lists! oooo
setwd("~/Projects/BcAt_RNAGWAS/data")
HiHerit <- read.csv("Wei_Bc/Wei2018_SDS5_TopHerit.csv")
LesCor <- read.csv("Wei_Bc/Wei2018_SDS6_CorLesBcTranscripts.csv")
BcNets <- read.csv("Wei_Bc/Wei2018_SDS7_BcNetworkGeneList.csv")
MixNets <- read.csv("Wei_Bc/Wei2018_SDS8_AtBcNetworkGeneList.csv")

require(reshape)
#HiHerit
names(HiHerit)[1] <- "GeneID"
HiHerit <- transform(HiHerit, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(HiHerit)[6] <- "Gene"
full.gene.hiherit <- merge(fg.sum, HiHerit, by = "Gene")
write.csv(full.gene.hiherit, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/SDS5_Bc_topGene.csv")
#LesCor
LesCor <- transform(LesCor, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(LesCor)[12] <- "Gene"
LesCor$Sig <- ifelse(LesCor$P.Value_All < 0.05, "Y", 
                     ifelse(LesCor$P.Value_Col0 < 0.05, "Y",
                            ifelse(LesCor$P.Value_coi1 < 0.05, "Y",
                                   ifelse(LesCor$P.Value_npr1 < 0.05, "Y", "N"))))
LesCor <- LesCor[LesCor$Sig=="Y",]
full.gene.lescor <- merge(fg.sum, LesCor, by = "Gene")
write.csv(full.gene.lescor, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/SDS6_Bc_topGene.csv")
#BcNets
BcNets <- transform(BcNets, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(BcNets)[13] <- "Gene"
full.gene.bcnets <- merge(fg.sum, BcNets, by = "Gene")
write.csv(full.gene.bcnets, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/SDS7_Bc_topGene.csv")
#MixNets
MixNets <- transform(MixNets, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(MixNets)[10] <- "Gene"
full.gene.mixnets <- merge(fg.sum, MixNets, by = "Gene")
write.csv(full.gene.mixnets, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/top10SNP/SDS8_Bc_topGene.csv")
#------------------------------------------------------------------------------------------------
