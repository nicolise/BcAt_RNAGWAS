#Nicole E Soltis
#split transcript data:
  #with significant Isolate fx
  #without sig Iso fx
  #plant?
#------------------------------------------------------------------------------------

#wei data
setwd("~/PhD/Literature/Collaborators/2018 - TBD - Botrytis transcriptional variation in response to Arabidopsis")
bcherit <- read.csv("SupplDS3_Pvalue.csv")
#keep necessary columns
bcherit <- bcherit[,c(1,21,26,31,37:39)]
