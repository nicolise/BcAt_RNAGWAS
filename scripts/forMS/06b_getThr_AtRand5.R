#Nicole E Soltis
#02/15/19

#-----------------------------------------------------------------------------
rm(list=ls())
#plot/ table: How often is max real p < max permut p?

#work on thresholding -- summarize across 5 permutations
#collect max value across 5 permuts for each SNP

setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At/")
#for now, don't care about annotating which phenotype was which (gene names for transcripts)
#do this if I only want ONE top SNP per transcript. But really, I want the top p value that any SNP ever gets, across any of my random transcripts.
mydat_r1 <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_RAND1_top1SNPsample.txt")
mydat_r1 <- mydat_r1[,-c(2)]
mydat_r1$randrun <- 1
mydatall <- mydat_r1
combdat <- mydat_r1[1,]
for (i in c(2:5)){
  mydat01 <- read.csv(paste("06_GEMMAsumm_RAND/col0_GEMMA_RAND",i,"_top1SNPsample.txt", sep=""))
  mydat01 <- mydat01[,-c(2)]
  mydat01$randrun <- i
  for (j in c(1:9267)){
    if(mydat01[j,12] < mydatall[j,12]) {
      combdat[j,] <- mydat01[j,]
    } else { combdat[j,] <- mydatall[j,]}
  }
  mydatall <- combdat
}
#min p value across all:
min(combdat$p_score) # 2.926053e-08
-log10(2.926053e-08)
write.csv(combdat, "06_GEMMAsumm_RAND/col0_GEMMA_maxRAND_1SNP.csv")
combdat <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_maxRAND_1SNP.csv")
#--------------------------------------------------------------
#plot max rand!
#Make plotting variables for snp
mydat <- combdat
mydat_plot <- mydat[order(mydat$chr, mydat$ps),]
mydat_plot$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr==i, ]$Index=mydat_plot[mydat_plot$chr==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr==i, ]$Index=mydat_plot[mydat_plot$chr==i, ]$ps+lastbase
  }
}

#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
#myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
myColors <- rep(c("darkgreen", "palegreen3"),9)
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/5xRand/AtCol0_top1SNP_MaxRand.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()

#---------------------------------------------------------------------------
#name transcripts for rand
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At/")
myPhenos <- NULL
nameddat <- NA
mydat <- NA
##check correct phenotype
myPhenos <- read.table("col0/02_GEMMA/binMAF20NA10.fam")
nameslist <- myPhenos[1,6:length(myPhenos)]
#check that V6 is a real phenotype
nameslist[1,1:10] #yes
myGenes <- as.data.frame(t(nameslist))
names(myGenes)[1] <- "Gene"
myGenes$pheno <- 1:nrow(myGenes)
#start with 1 file. Loop over all files in directory later
##change these out to annotate all SNP summaries per GEMMA run
mydat <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_maxRAND_1SNP.csv")
#, row.names=NULL
names(mydat)
nameddat <- merge(mydat, myGenes, by = "pheno")
##match name here
write.csv(nameddat, "05_GEMMAsumm/GeneNames/col0_GEMMA_top100SNPsample.csv")

#--------------------------------------------------------------------------

#plot/ table: How often is max real p < max permut p?
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At/")
randdat <- read.csv("06_GEMMAsumm_RAND/col0_GEMMA_maxRAND_top1SNPsample_named.csv")
#need annotated gene names ("phenotype") to match random to real data

mydat01 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")

randmerge <- randdat[,c("Gene","chr","ps", "p_score","randrun")]
names(randmerge) <- c("Gene","rand_chr","rand_pos", "rand_p","rand_run")
hotspt <- merge(mydat01, randmerge, by="Gene")
hotspt$DvR <- hotspt$rand_p - hotspt$p_score
hist(hotspt$DvR)
hotspt$mygroup <- ifelse(hotspt$DvR > 0, "NonSig", "Sig")
#add indexing now
mydat_plot <- hotspt
mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
mydat_plot$Index.s = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr==i, ]$Index.s=mydat_plot[mydat_plot$chr==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr==i, ]$Index.s=mydat_plot[mydat_plot$chr==i, ]$ps+lastbase
  }
}

#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
mydat_plot$chr.sig <- as.factor(paste(mydat_plot$chr, mydat_plot$mygroup, sep="."))
levels(mydat_plot$chr.sig)
myColors <- c("navyblue","darkred", "royalblue1","indianred1",  "navyblue","darkred",  "royalblue1","indianred1", "navyblue","darkred", "royalblue1", "indianred1", "navyblue","darkred", "royalblue1","indianred1",  "navyblue","darkred",  "royalblue1","indianred1",  "royalblue1","indianred1", "navyblue","darkred",  "royalblue1","indianred1", "navyblue","darkred", "royalblue1","indianred1","navyblue", "darkred",  "royalblue1","indianred1", "navyblue","darkred",  "royalblue1", "indianred1")
names(myColors) <- levels(mydat_plot$chr.sig)
colScale <- scale_colour_manual(name = "Sig",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_top1SNP_realOverRand.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index.s, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(mydat_plot$chr.sig),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()

#summarize real SNP vs. permuted SNP: what % of the time is real data more sig?
table(mydat_plot$mygroup)
6405/(2862+6405) #69%
#------------------------------------------------------------------------------------
#additional bits using simulation threshold
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachBc_At")
mydat01 <- read.csv("05_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")
#Make plotting variables for snp
mydat <- mydat01
mydat_plot <- mydat[order(mydat$chr, mydat$ps),]
mydat_plot$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr==i, ]$Index=mydat_plot[mydat_plot$chr==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr==i, ]$Index=mydat_plot[mydat_plot$chr==i, ]$ps+lastbase
  }
}

#summarize within each SNP - # of transcript hits
mydat_summ <- mydat_plot[,c("chr","ps","p_score","Gene","Index")]
mydat_summ_ngene <- aggregate(Gene ~ Index, data = mydat_summ, FUN = function(x){NROW(x)})
#now add SNP data back on, matching by Index.s
mydat_summ_ngene <- merge(mydat_summ_ngene, mydat_summ[,c("chr","ps","Index")], by="Index")
#remove duplicate rows
mydat_summ_ngene <- unique(mydat_summ_ngene)
mydat_plot <- mydat_summ_ngene

library(ggplot2)
#create a custom color scale
myColors <- rep(c("darkgreen", "palegreen3"), 9)
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chr",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
## check name depending on df
jpeg("plots/Manhattans/AtCol0_top1SNP_GeneCounts.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=Gene))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(mydat_plot$chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="Number of Genes")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()

#here, with dashed line for simulation threshold at 5 genes overlap
setwd("~/Projects/BcAt_RNAGWAS")
## check name depending on df
jpeg("plots/Manhattans/AtCol0_top1SNP_GeneCounts_5thr.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=Gene))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(mydat_plot$chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="Number of Genes")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    geom_hline(aes(yintercept=5), linetype=2)+
    expand_limits(y=0)
)
dev.off()
#----------------------------------------------------------------
setwd("~/Projects/BcAt_RNAGWAS")
mydat_c9 <- mydat[mydat$chr==9,]
jpeg("plots/Manhattans/5xRand/AtCol0_chr9_top1SNP.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_c9, aes(x=ps, y=(-log10(p_score))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Distance (Mb)", breaks=c(0,1e+06,2e+06,3e+06),labels=c(0,1,2,3))+
    expand_limits(y=0)
)
dev.off()