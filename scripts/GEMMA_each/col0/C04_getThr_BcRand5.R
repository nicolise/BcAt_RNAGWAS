#Nicole E Soltis
#11/05/18

#-----------------------------------------------------------------------------
#work on thresholding -- summarize across 5 permutations
#collect max value across 5 permuts for each SNP

rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
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
min(combdat$p_score) #2.949698e-08 
-log10(2.949698e-08 )
write.csv(combdat, "06_GEMMAsumm_RAND/col0_GEMMA_maxRAND_1SNP.csv")

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
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
myColors <- c("navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1")
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/5xRand/BcCol0_top1SNP_MaxRand.jpg", width=8, height=5, units='in', res=600)
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

#--------------------------------------------------------------------
#mini manhattan: focus in on chr 1
#plot max rand on Chr 1 only. Label Mb X axis. 
setwd("~/Projects/BcAt_RNAGWAS")
mydat_c1 <- mydat_plot[mydat_plot$chr==1,]
jpeg("plots/Manhattans/5xRand/BcCol0_chr1_top1SNP_maxrand.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_c1, aes(x=Index, y=(-log10(mydat_c1$p_score))))+
    theme_bw()+
    colScale+
    #used stroke = 0 for top 10, not top 1
    #, stroke=0
    geom_point(aes(color = factor(chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Distance (Mb)", breaks=c(0,1e+06,2e+06,3e+06,4e+06),labels=c(0,1,2,3,4))+
    expand_limits(y=0)
)
dev.off()

#---------------------------------------------------------------------------

#plot/ table: How often is max real p < max permut p?
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
randdat <- read.csv("06_GEMMAsumm_RAND/TranscNames/col0_GEMMA_top1SNPsample_genes.csv")
#need annotated gene names ("phenotype") to match random to real data

mydat01 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")

randmerge <- randdat[,c("Gene","chr","ps", "p_score","randrun")]
names(randmerge) <- c("Gene","rand_chr","rand_pos", "rand_p","rand_run")
hotspt <- merge(mydat01, randmerge, by="Gene")
hotspt$DvR <- hotspt$rand_p - hotspt$p_score
hist(hotspt$DvR)
hotspt$mycolor <- ifelse(hotspt$DvR > 0, "NonSig", "Sig")

mydat_plot <- hotspt
#plot by SNP location! (disregard transcript location- just want hotspots)
library(ggplot2)
#create a custom color scale
myColors <- c("navyblue", "royalblue1")
names(myColors) <- levels(mydat_plot$mycolor)
colScale <- scale_colour_manual(name = "Sig",values = myColors)

mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
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

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_top1SNP_realOverRand.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=(-log10(mydat_plot$p_score))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(mydat_plot$mycolor),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="-log10(p)", breaks=c(0,2.5,5,7.5), labels=c("0","2.5","5.0","7.5"), limits=c(0,8.75))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()

#could also plot this as distance -- set all positive values to zero, then plot -1*x
mydat_plot$DvR <- ifelse(mydat_plot$DvR > 0, 0, mydat_plot$DvR)

library(ggplot2)
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
myColors <- c("navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1","navyblue", "royalblue1")
names(myColors) <- levels(mydat_plot$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/Manhattans/BcCol0_top1SNP_realOverRand_dist.jpg", width=8, height=5, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=((-1*DvR))))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(mydat_plot$chr),alpha=0.001))+
    labs(list( title=NULL))+
    theme(legend.position="none")+
    scale_y_continuous(name="Distance Observed vs. Random")+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))+
    expand_limits(y=0)
)
dev.off()
#---------------------------------------------------------------------
#other version: pair based on SNP location, not based on transcript location
#currently, can assume that most TOP SNPs differ between each rand list. So just combine all (rbind) and then for duplicated values of chr.snp, take the top (lowest p) value


