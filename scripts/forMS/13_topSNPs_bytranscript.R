#Nicole E Soltis
#03/01/19

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#Bc transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/04_GEMMAout/
#kmat script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A03_GEMMA_kmat.sh
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc/A04_runGEMMA_kmat_3genos.sh

#Need the following:
#median/ distribution of number of sig SNPs per transcript
#look for extremely low p values = evidence for cis effects

rm(list=ls())
#setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")
setwd("D:/BcAt_RNAGWAS")

#get list of phenos with many SNP > Thr
mysnplistB <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/Bc_phenos_manySNPovrThr.csv")
mysnplistA <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/At_phenos_manySNPovrThr.csv")
totnumsnpB <- read.csv("GEMMA_eachAt_Bc/07_TopSNPs/Bc_numphenosoverThr.csv")
#-----------------------------------------------------------------------------------
#calculate median SNP > thr 

#-----------------------------------------------------------------------------------
#look for extreme p-values as evidence for cis
rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/")
top100B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top100SNPsample.txt")
top1B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top1SNPsample.txt")

setwd("~/Projects/BcAt_RNAGWAS/data")
top100B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top100SNPsample.txt")
top1B <- read.csv("GEMMA_eachAt_Bc/06_GEMMAsumm/col0_GEMMA_top1SNPsample.txt")
top100A <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top100SNPsample_genes.csv")
top1A <- read.csv("GEMMA_eachBc_At/05_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top1SNPsample_genes.csv")

#top 100 Bc
hist(-log10(top100B$p_score))
min(top100B$p_score)
hi100B <- top100B[top100B$p_score < 1e-6,]
hist(hi100B$p_score)
hist(-log10(hi100B$p_score))

#top 1 Bc
hist(-log10(top1B$p_score))
min(top1B$p_score)
hi1B <- top1B[top1B$p_score < 1e-6,]
hist(hi1B$p_score)
hist(-log10(hi1B$p_score))
quantile(hi1B$p_score, c(0.00, 0.01, 0.05, 0.25, 0.5))

#top 100 At
hist(-log10(top100A$p_score))
min(top100A$p_score)
hi100A <- top100A[top100A$p_score < 1e-6,]
hist(hi100A$p_score)
hist(-log10(hi100A$p_score))

#top 1 Bc
hist(-log10(top1A$p_score))
min(top1A$p_score)
hi1A <- top1A[top1B$p_score < 1e-6,]
hist(hi1A$p_score)
hist(-log10(hi1A$p_score))
quantile(hi1A$p_score, c(0.00, 0.01, 0.05, 0.25, 0.5))

library(ggplot2)
ggplot(hi1A, aes(y=-log10(p_score), x=1)) + 
  geom_violin()

ggplot(hi1B, aes(y=-log10(p_score), x=1)) + 
  geom_violin()

boxplot(-log10(hi1A$p_score), ylim=c(0, 9))
boxplot(-log10(hi1B$p_score), ylim=c(0, 9))

boxplot(-log10(hi100A$p_score), ylim=c(0, 9))
boxplot(-log10(hi100B$p_score), ylim=c(0, 9))
#use the Tukey’s method to identify the outliers ranged above and below the 1.5*IQR
boxplot.stats(-log10(hi1B$p_score))$out
boxplot.stats(-log10(hi1A$p_score))$out

#better plots - histograms
library(ggplot2)
ggplot(hi1B, aes(x=-log10(p_score))) + geom_histogram(color="navyblue", fill="royalblue1") + 
  geom_vline(aes(xintercept=median(-log10(hi1B$p_score))), color="navyblue", linetype="dashed") + 
  theme_bw()

ggplot(hi1A, aes(x=-log10(p_score))) + geom_histogram(color="darkgreen", fill="palegreen3") + 
  geom_vline(aes(xintercept=median(-log10(hi1A$p_score))), color="darkgreen", linetype="dashed") + 
  theme_bw()

#check whether the top 50 p-value SNPs are in cis
#or top 10, but wahtever
top1B.max.p <- top1B[order(top1B$p_score),]
top1B.max.p <- top1B.max.p[1:50,]
top1B.max.eff <- top1B[order(top1B$beta),]
top1B.max.eff <- top1B.max.eff[1:50,]

write.csv(top1B.max.p, "GEMMA_eachAt_Bc/07_TopSNPs/Bc_topSNPs_p_cis.csv")
write.csv(top1B.max.eff, "GEMMA_eachAt_Bc/07_TopSNPs/Bc_topSNPs_beta_cis.csv")
#-------------------------------------------------------------------------------
#now annotate these phenotypes to genes to tell if in cis or not
top1B.max.p$group <- "max.p"
top1B.max.eff$group <- "max.eff"
top1B <- rbind(top1B.max.p, top1B.max.eff)
#first name phenotypes with genes
#setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachAt_Bc")
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/02_GEMMA")
myphenolookup <- read.csv("PhenoToGene.csv")
myphenolookup <- myphenolookup[,-c(1)]
top1B.gen <- merge(top1B, myphenolookup, by="pheno")

#-------------------------------------------------------------------
mydat <- top1B.gen

setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
num.genes <- my.gtf[unique(my.gtf$V12),]
my.gtf <- my.gtf[,1:14]
my.gtf$midgene <- (my.gtf$V4 + my.gtf$V5)/2

mydat.genes <- mydat
mydat.genes$GeneNoTranscript <- gsub("\\.[0-9]$", '', mydat.genes$Gene)
my.gtf.genes <- my.gtf[,c("V1","V3","V4","V5","V10","V12","V14")]
my.gtf.genes$V12 <- gsub("^gene:","", my.gtf.genes$V12)
my.gtf.genes$V10 <- gsub("^transcript:","", my.gtf.genes$V10)
names(my.gtf.genes) <- c("chr", "exon","start","stop","transcript","Gene","GeneName")
my.gtf.genes$chr <- gsub("^Chromosome","", my.gtf.genes$chr)
names(mydat.genes)[17] <- "transcript"
names(mydat.genes)[18] <- "Gene"

library(dplyr)
my.gtf_genesumm <- my.gtf.genes %>%
  group_by(transcript, GeneName, chr) %>%
  summarize(tstart = min(start, na.rm=TRUE),
            tstop = max(stop, na.rm=TRUE))
my.gtf_genesumm <- as.data.frame(my.gtf_genesumm)
my.gtf_genesumm$tmid <- (my.gtf_genesumm$tstop + my.gtf_genesumm$tstart)/2
names(my.gtf_genesumm)[3] <- "chr.t"

mydat_cis <- merge(mydat.genes, my.gtf_genesumm, by="transcript")

mydat_iscis <- mydat_cis[mydat_cis$chr==mydat_cis$chr.t,]

mydat_iscis$ps
mydat_iscis$tmid

mydat_cis$dist <- ifelse(mydat_cis$chr!=mydat_cis$chr.t, "trans", abs(mydat_cis$ps - mydat_cis$tmid))

#copied to Supplemental Table 2
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc")
write.csv(mydat_cis, "07_TopSNPs/Bc_topSNPscis.csv")
