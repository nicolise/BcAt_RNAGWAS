#Nicole E Soltis
#06/05/18
#Focus on phenotypes in networks of interest
#----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
setwd("~/Projects/BcAt_RNAGWAS/")
Phenos <- read.csv("data/allreads_bigRR/BO5.10/01_prepFiles/lsmeans_zscale_allreads.csv")
Phenos <- Phenos[,-c(1)]

#members of network 5: 
MyNet5List <- c("Bcin01g11450.1", "Bcin01g11460.1", "Bcin01g11470.1", "Bcin01g11470.2", "Bcin01g11480.1", "Bcin01g11490.1", "Bcin01g11500.1", "Bcin01g11520.1", "Bcin01g11530.1", "Bcin01g11550.1")

#bot network
#from Wei Supp_DS7
MyBotList <- c("Bcin12g06370.1", "Bcin12g06380.1", "Bcin12g06390.1", "Bcin12g06400.1", "Bcin12g06410.1", "Bcin12g06420.1", "Bcin12g06430.1")

#boa network
#from Wei Supp_DS7
MyBoaList <- c("Bcin01g00010.1", "Bcin01g00020.1", "Bcin01g00030.1", "Bcin01g00040.1", "Bcin01g00050.1", "Bcin01g00060.1", "Bcin01g00070.1", "Bcin01g00080.1", "Bcin01g00090.1", "Bcin01g00100.1", "Bcin01g00110.1", "Bcin01g00120.1", "Bcin01g00130.1")

MyNetsList <- c(MyBotList, MyBoaList, MyNet5List)

#make a dataframe of members of each network
MyNets <- as.data.frame(MyNetsList)
names(MyNets)[1] <- "Gene"
MyNets$Cluster <- ifelse(MyNets$Gene %in% MyBoaList, "BOA", ifelse(MyNets$Gene %in% MyBotList, "BOT", "NET5"))
write.csv(MyNets, "data/BcBotGWAS/02_MatchGenos/NetworkMembers.csv")

PhenosNet <- Phenos[,names(Phenos) %in% c("Isolate",MyNetsList)]

write.csv(PhenosNet, "data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")
