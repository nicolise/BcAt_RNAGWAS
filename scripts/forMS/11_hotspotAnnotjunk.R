setwd("~/Projects/BcAt_RNAGWAS/paper/Tables/formatted")
mygenes <- read.csv("extra/HotspotTargets_At_allgenes.csv")
mygo <- read.csv("extra/HotspotTargets_At_GOlist.csv")
names(mygenes)[1] <- "Gene"
names(mygo)[1] <- "Gene"
mygo <- merge(mygo, mygenes, by="Gene")
write.csv(mygo, "extra/HotspotTargets_At_GOlist_peaksTagged.csv")

commongo <- as.data.frame(table(mygo$Function))
names(commongo)[1] <- "fulldata"
blah <- split(mygo, mygo$eQTLpeak)
for (i in 1:20){
  bleh <- as.data.frame(blah[i])
  names(bleh)[2] <- "Function"
  bleh$Function <- droplevels(bleh$Function)
  bleh <- as.data.frame(table(bleh$Function))
  
}