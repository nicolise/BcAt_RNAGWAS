#Nicole E Soltis
#06/08/18
#03_bigRR_cistrans
#--------------------------------------------------------------------
#extract BotBoaNet5 bigRR data to subfolder
# try post-hoc meta-analysis across phenotypes
#check cis vs. trans SNP effect estimates
#later repeat this for GEMMA
rm(list = ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
#here are the original phenotypes for these SNPs
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")

#now extract relevant GWAS data from here...
#from data/allreads_bigRR/B05.10/03_bigRRout_partial (if missing check .tar.gz)
#to data/allreads_bigRR/B05.10/04_NetSubset