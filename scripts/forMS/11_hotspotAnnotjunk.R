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
  summ <- as.data.frame(table(bleh[,2]))
  commongo <- cbind(commongo, summ[,2])
  names(commongo)[2+i] <- paste0(names(bleh)[2])
}
write.csv(commongo, "extra/HotspotTargets_At_GOlist_commonterms.csv")
