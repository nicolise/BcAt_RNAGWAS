#Nicole E Soltis
#01/11/19
#grab genes under hotspots, look for network overlap/ functions

#---------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS")
TopSNP.genes <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/BcAt_topSNP_Genes_ed.csv")
SigHots <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/BcAt_sigGenes.csv")
SigHot.genes <- TopSNP.genes[TopSNP.genes$Gene %in% SigHots$Gene,]

#get Bc transcripts under sig hotspots
BcDat <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")

#get At transcripts under sig hotspots
AtDat <- read.csv("data/GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top1SNPsample_genes.csv")

BcDat$chr_snp <- paste(BcDat$chr, BcDat$ps, sep="_")
AtDat$chr_snp <- paste(AtDat$chr, AtDat$ps, sep="_")

#for both of these, "Gene" means Transcript
BcDat_sig <- BcDat[BcDat$chr_snp %in% SigHot.genes$chr_snp,]
AtDat_sig <- AtDat[AtDat$chr_snp %in% SigHot.genes$chr_snp,]

SH.genes <- SigHot.genes[,c("Gene","chr_snp")]
SH.genes <- unique(SH.genes) #cool, 26 hotspots. manageable.
names(SH.genes)[1] <- "HotSpotNearestGene"

BcDat_sig <- merge(BcDat_sig, SH.genes, by="chr_snp")
AtDat_sig <- merge(AtDat_sig, SH.genes, by="chr_snp")

BcDat_summ <- BcDat_sig[,c("Gene","HotSpotNearestGene")]
AtDat_summ <- AtDat_sig[,c("Gene","HotSpotNearestGene")]
#write out dataframes
setwd("~/Projects/BcAt_RNAGWAS")
write.csv(BcDat_summ, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_genesUnderSigHotspots.csv")
write.csv(AtDat_summ, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_genesUnderSigHotspots.csv")

#-----------------------------------------------------------------------------------------------
#look for function annotations for the Botrytis genes!
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")

BcDat_summ <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_genesUnderSigHotspots.csv")
AtDat_summ <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_genesUnderSigHotspots.csv")
BcDat_summ <- BcDat_summ[,-c(1)]

funclist <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_HotSpot_Botportal.csv")
names(funclist)[1] <- "Gene"

Bc_fn <- merge(BcDat_summ, funclist, by="Gene", all.x=TRUE)
write.csv(Bc_fn,"data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_Hotspot_funcannot.csv")

#and try connecting this list to Vivian's Bc lists! oooo
setwd("~/Projects/BcAt_RNAGWAS/data")
HiHerit <- read.csv("Vivian_Bc/Wei2018_SDS5_TopHerit.csv")
LesCor <- read.csv("Vivian_Bc/Wei2018_SDS6_CorLesBcTranscripts.csv")
BcNets <- read.csv("Vivian_Bc/Wei2018_SDS7_BcNetworkGeneList.csv")
MixNets <- read.csv("Vivian_Bc/Wei2018_SDS8_AtBcNetworkGeneList.csv")


require(reshape)
#HiHerit
names(HiHerit)[1] <- "GeneID"
HiHerit <- transform(HiHerit, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(HiHerit)[6] <- "Gene"
Bc.gene.hiherit <- merge(BcDat_summ, HiHerit, by = "Gene")
write.csv(Bc.gene.hiherit, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_SDS5_HotSpots.csv")
#LesCor
LesCor <- transform(LesCor, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(LesCor)[12] <- "Gene"
LesCor$Sig <- ifelse(LesCor$P.Value_All < 0.05, "Y", 
                     ifelse(LesCor$P.Value_Col0 < 0.05, "Y",
                            ifelse(LesCor$P.Value_coi1 < 0.05, "Y",
                                   ifelse(LesCor$P.Value_npr1 < 0.05, "Y", "N"))))
LesCor <- LesCor[LesCor$Sig=="Y",]
Bc.gene.lescor <- merge(BcDat_summ, LesCor, by = "Gene")
write.csv(Bc.gene.lescor, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_SDS6_HotSpots.csv")
#BcNets
BcNets <- transform(BcNets, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(BcNets)[13] <- "Gene"
Bc.gene.bcnets <- merge(BcDat_summ, BcNets, by = "Gene")
write.csv(Bc.gene.bcnets, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_SDS7_HotSpots.csv")
#MixNets
MixNets <- transform(MixNets, a = colsplit(GeneID, split = "\\.", names = c('Gene', 'Transc')))
names(MixNets)[10] <- "Gene"
Bc.gene.mixnets <- merge(BcDat_summ, MixNets, by = "Gene")
write.csv(Bc.gene.mixnets, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_SDS8_HotSpots.csv")

#for At
#MixNets
At.gene.mixnets <- merge(AtDat_summ, MixNets, by = "Gene")
write.csv(At.gene.mixnets, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_SDS8_HotSpots.csv")

#and her At genes list for networks
setwd("~/Projects/BcAt_RNAGWAS/data")
#figure 5: gene co-expression networks
#supp fig 5: gene co-expression archictecture Col-0
#table 5: hub genes and bottlenecks
AtHub <- read.csv("Vivian_At/Table5_HubGenes.csv")
#SDS 6: genes condensed in co-expression networks
AtGCN <- read.csv("Vivian_At/SDS6_GCN.csv")
#SDS 8: tables of top 5% genes from PCA
AtPCA <- read.csv("Vivian_At/SDS8_top5pctGenes.csv")

#remove empty column
AtDat_summ <- AtDat_summ[,-c(1)]
#AtHub
names(AtHub)[1:2] <- c("Network","Gene")
At.gene.hub <- merge(AtDat_summ, AtHub, by="Gene") #50 genes
write.csv(At.gene.hub, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_2016_Table5_Hotspots.csv")
#AtGCN
names(AtGCN)[1:2] <- c("NetworkMember","Gene")
At.gene.gcn <- merge(AtDat_summ, AtGCN, by="Gene") #529 genes
write.csv(At.gene.gcn, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_2016_SDS6_Hotspots.csv")
#just gene hits
At.gene.gcn.sm <- At.gene.gcn[,c("Gene","HotSpotNearestGene","Gene.Name")]
At.gene.gcn.sm <- unique(At.gene.gcn.sm)
write.csv(At.gene.gcn.sm, "GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_2016_SDS6_Hotspots_summary.csv")
#AtPCA
names(AtPCA)[1] <- "Gene"
At.gene.pca <- merge(AtDat_summ, AtPCA, by="Gene")#233 genes
write.csv(At.gene.pca, "GEMMA_eachAt_Bc/07_TopSNPS/BcAt_permut/At_2016_SDS8_Hotspots.csv")
#just gene hit counts
At.gene.pca.sm <- At.gene.pca[,c("Gene","HotSpotNearestGene","Gene.Name")]
At.gene.pca.sm <- unique(At.gene.pca.sm) #146
write.csv(At.gene.pca.sm, "GEMMA_eachAt_Bc/07_TopSNPS/BcAt_permut/At_2016_SDS8_Hotspots_summary.csv")
#------------------------------------------------------------------------------------------------
#now, split out into lists according to which hotspot each transcript is linked to
#first, remove extra levels of factor
BcDat_sig$HotSpotNearestGene <- droplevels(BcDat_sig$HotSpotNearestGene)
BcSplit <- split( BcDat_sig , f = BcDat_sig$HotSpotNearestGene)

AtDat_sig$HotSpotNearestGene <- droplevels(AtDat_sig$HotSpotNearestGene)
AtSplit <- split( AtDat_sig, f = AtDat_sig$HotSpotNearestGene)

