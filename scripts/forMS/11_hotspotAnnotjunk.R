setwd("~/Projects/BcAt_RNAGWAS/paper/Tables/formatted")
mygenes <- read.csv("HotspotTargets_At_allgenes.csv")
mygo <- read.csv("extra/HotspotTargets_At_GOlist.csv")
names(mygenes)[1] <- "Gene"
names(mygo)[1] <- "Gene"
mygo <- merge(mygo, mygenes, by="Gene")
write.csv(mygo, "HotspotTargets_At_GOlist_peaksTagged.csv")
