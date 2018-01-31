#reformat bigRR output data
#Nicole E Soltis

#--------------------------------------------------------
rm(list=ls())
library(tidyr)
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcSolGWAS/data/GWAS_files/04_bigRRoutput/trueMAF20_20NA/")
#Import data
HEM.plotdata <- read.csv( "SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv") 
HEM.thresh <- read.csv("SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")
#read in to 06_plots