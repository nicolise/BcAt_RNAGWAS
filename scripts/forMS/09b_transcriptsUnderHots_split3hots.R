#Nicole E Soltis
#01/11/19
#grab genes under hotspots, look for network overlap/ functions
#split 3 hotspots by SNP, not gene level: Bcin09g06590 and Bcin12g00330 (Bc), Bcin12g02130 (At)
#Bcin09g06590 is 9_2330312 and 9_2334368 - 4kb apart, may be unique loci
#Bcin12g00330 is 12_115491 and 12_115511 - 20 bp apart, almost certainly splitting 1 locus
#Bcin12g02130 is 12_758420 and 12_760499 - 2kb apart, may be unique loci

#---------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS")
TopSNP.genes <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/BcAt_topSNP_Genes_ed.csv")
SigHots <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/BcAt_sigGenes.csv")
SigHot.genes <- TopSNP.genes[TopSNP.genes$Gene %in% SigHots$Gene,]

myhots_list <- c("9_2330312","9_2334368", "12_115491", "12_115511", "12_758420", "12_760499")

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

BcDat_summ <- BcDat_sig[,c("Gene","chr_snp")]
AtDat_summ <- AtDat_sig[,c("Gene","chr_snp")]

#write out dataframes
setwd("~/Projects/BcAt_RNAGWAS")
write.csv(BcDat_summ, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_genesUnderSigHotspots_SNPlvl.csv")
write.csv(AtDat_summ, "data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/At_genesUnderSigHotspots_SNPlvlv.csv")
BcDat_summ <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/Bc_genesUnderSigHotspots_SNPlvl.csv")

#modify existing DF to add in chr_snp information on hotSNPs
setwd("~/Projects/BcAt_RNAGWAS")
myAta <- read.csv("paper/Tables/formatted/extra/Table N3a. HotspotTargets_At_annotation.csv")
myAtb <- read.csv("paper/Tables/formatted/extra/Table N3b. HotspotTargets_At_annotation.csv")
myBc <- read.csv("paper/Tables/formatted/extra/Table N2. HotspotTargets_Bc_annotation.csv")
myBc <- merge(myBc, BcDat_summ, by="Gene")
write.csv(myBc, "paper/Tables/formatted/extra/Table N2. HotspotTargets_Bc_annotation_ed.csv")
names(myAta)[1] <- "Gene"
myAta <- merge(myAta, AtDat_summ, by="Gene")
write.csv(myAta, "paper/Tables/formatted/extra/Table N3a. HotspotTargets_At_annotation_ed.csv")
names(myAtb)[1] <- "Gene"
myAtb <- merge(myAtb, AtDat_summ)
write.csv(myAtb, "paper/Tables/formatted/extra/Table N3b. HotspotTargets_At_annotation_ed.csv")
