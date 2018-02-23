#Nicole E Soltis
#02/21/18
#02_FAMaddPhenos

#----------------------------------------------------------------------------
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
Phenos <- read.csv("data/allreadsGWAS/BO5.10/01_prepFiles/lsmeans_zscale_allreads.csv")
Phenos <- Phenos[,-c(1)]

#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("data/B05_GEMMA/01_PLINK/dp_charMAF20_10NA.fam")