#Nicole E Soltis
#11/28/18
#linear models for each BOA gene to test effect of local genetic variation on gene expression
#------------------------------------------------------------------------------
rm(list=ls())
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_boa_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_boa_fix_named.csv", na.strings=c("","NA"))
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/")
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
PhenosNet <- PhenosNet[,-c(1)]
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")
myclusters <- read.csv("data/BcBotGWAS/BotBoaNet5_pvclust.csv")
names(myclusters)[1] <- "Isolate"
#get just BOA gene regions for plot part C
#end line 158
rm(list=ls())
#gtf location: 
setwd("~/Projects/BcGenome")
my.gtf <- read.table("data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE)
my.gtf <- my.gtf[,c(1:8,10,12,14)]
my.gtf.c1 <- my.gtf[my.gtf$V1=="Chromosome1",]
my.gtf.boa <- my.gtf[1:158,]

setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_boa_named <- read.csv("02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa_rmcalls.csv", na.strings=c("","NA"))
#try SNP list without rmcalls to see why there are gaps between genes
mySNP_boa_b <- read.csv("02b_Haploview/BOA_deletion/binMAF20NA10_chr1_boa_wgenes.csv")
#gaps (with no SNPs) between genes are due to lack of SNP calls, not our dropped SNPs
mySNP_boa <- mySNP_boa_named
snplist <- as.numeric(as.character(t(as.vector(mySNP_boa[2,2:length(mySNP_boa)]))))

#now, mean phenotype plots by cluster membership
#assign isolates to clusters
BoaNets <- MyNets[MyNets$Cluster=="BOA",]
BoaPhenos <- PhenosNet[,names(PhenosNet) %in% BoaNets$Gene]
BoaPhenos <- cbind(PhenosNet$Isolate, BoaPhenos)
names(BoaPhenos)[1] <- "Isolate"
BoaPhenos$mean.Pheno <- rowMeans(BoaPhenos[,c(2:length(BoaPhenos))])
#need a new file of cluster membership
BoaClusts <- myclusters[,c("Isolate","boa_pv_gene_indel")] #can also do with boa_pv_gene
BoaPhenosPV <- merge(BoaPhenos,BoaClusts, by = "Isolate")

#linear model of gene expression 
mySNP <- mySNP_boa_named[2:98,c(1,seq(2, 469, by=2))]
names(mySNP) <- as.character(unlist(mySNP[1,]))
names(mySNP)[1] <- "Isolate"
mySNP <- mySNP[2:97,]
#only keep every other column, plus column 1
myphenos <- BoaPhenosPV

mymoddat <- merge(myphenos, mySNP, by="Isolate")

mymodouts <- as.data.frame(NULL)
#trying SNPs as fixed effects, not sure
#16 is boa_pv_gene_indel, which accounts for the isolate groupings including the deletion.
for (i in c(2:14)){
mymodfactors <- mymoddat[,c(i,16:250)]
names(mymodfactors)[1] <- "pheno"
mymod <- lm(mymodfactors[,1] ~ . - pheno, data=mymodfactors)
blah <- as.data.frame(anova(mymod))
blah$fdr.p <- p.adjust(blah$`Pr(>F)`, method="fdr")
blah$Cats <- ifelse(blah$fdr.p > 0.05, "NS", 
                    ifelse(blah$fdr.p > 0.01, "<0.05",
                           ifelse(blah$fdr.p > 0.001, "<0.01","<0.001")))
if(i == 2) {mymodouts <- as.data.frame(blah$Cats)} else { mymodouts <- cbind(mymodouts, blah$Cats)}
names(mymodouts)[i-1] <- names(mymoddat[i])
}
rownames(mymodouts) <- rownames(blah)
