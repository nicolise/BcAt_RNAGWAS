#Nicole E Soltis
#05/10/18

#reminder of where I stored JAC metabolite data (GLS and camalexin)
#-----------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")

myGLS <- read.csv("Metabolites/GLSTotal_JAC.csv")
myCam <- read.csv("Metabolites/BcAtGWAS_TotalCam_JAC_avgCam.csv")
