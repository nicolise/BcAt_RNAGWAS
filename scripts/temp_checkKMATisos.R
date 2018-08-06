#check on kmat dimensions
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc")
myped <- read.table("01_PLINK/dpbinMAF20NA10.ped")
myfam <- read.table("01_PLINK/binMAF20NA10.fam") #only 95 isos
#mymap <- read.table("01_PLINK/dpbinMAF20NA10.map") #does not include isolate IDs
