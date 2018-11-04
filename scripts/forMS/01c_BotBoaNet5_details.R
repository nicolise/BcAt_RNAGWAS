#Nicole E Soltis
#get details about Bot, Boa, Net5 networks for text

#-----------------------------------------------------------------------------
rm(list=ls())
setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/")
#my nets gives gene members for each network
MyNets <- read.csv("data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")

setwd("C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm")
mySNP_bot_named <- read.csv("02b_Haploview/binMAF20NA10_chr12_bot_recrop_named.csv", na.strings=c("","NA"))
mySNP_boa_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_boa_fix_named.csv", na.strings=c("","NA"))
mySNP_net5_named <- read.csv("02b_Haploview/binMAF20NA10_chr1_net5_fix_named.csv", na.strings=c("","NA"))

#get start and end of each network... check mysnp_named files
#BOT = Bcin12g06370.1 to Bcin12g06430.1 
#start = 2217400, end = 2243357
2243357 - 2217400

#BOA = Bcin01g00030.1 to Bcin01g00130.1
#start = 5464, end = 61263
61263 - 5464
#BOA deletion boundaries?
#from excel chart: file:///C:/Users/nesol/Documents/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bclsm/02b_Haploview/BOA_deletion/TableX1_binMAF20NA10_chr1_boa_rmcalls.xlsx 
#start = 4029, end = 82614
82614 - 4029

#NET5 = Bcin01g11450.1 to Bcin01g11550.1
#start = 4026703, end = 4073215 
4073215 - 4026703
