#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")

#just Manhattan plots for individual plant genotypes
############################################################################
###Plotting the HEM results

#NEED TO CHECK CONTIGS FOR THIS

#Load plotting package
library(ggplot2)
library(grid)
library(dplyr)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata.og <- read.csv("data/BcBOT_MAF20_HEM.PlotFormat.csv")

HEM.plotdata <- HEM.plotdata.og

HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
HEM.thresh <- read.csv("data/BcBOT_MAF20_HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1)]


TH95pos <- HEM.thresh[1,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}
TH95neg <- HEM.thresh[5,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}
TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}
TH999pos <- HEM.thresh[4,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}
TH999neg <- HEM.thresh[8,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

#create a custom color scale
myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#get midpoint positions per chromosome
((max(HEM.plotdata[ which(HEM.plotdata$Chrom=='1'),]$Index) - min(HEM.plotdata[ which(HEM.plotdata$Chrom=='1'),]$Index))/2+min(HEM.plotdata[ which(HEM.plotdata$Chrom=='1'),]$Index))

#get length per chromosome segment
max(HEM.plotdata[which(HEM.plotdata$Chrom.Seg.F=='16.7'),]$Index) - min(HEM.plotdata[which(HEM.plotdata$Chrom.Seg.F=='16.7'),]$Index)

#this is the correct list for NA20
#c(1677875, )

#try: keep top SNPs only
for (i in c(4:24)){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}

#for top 1000 only
# for (i in c(4:24)){
#   assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
#   assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), tail(arrange(get(paste("HEMneg.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
# }

#combine pos and neg by group
for (i in c(4:15)){
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
}

#then combine
for (i in c(4:15)){
  mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[5] <- "Effect"
  assign(mydf, renamedf)
  myblob <- rep(names(HEM.plotdata[i]), nrow(get(mydf)))
  assign(mydf, cbind(get(mydf), Trait = myblob))
}

HEM.topSNPs <- rbind(HEM.Bcin12g06370.1.coi.1,HEM.Bcin12g06380.1.coi.1,HEM.Bcin12g06380.1.col.0,HEM.Bcin12g06390.1.coi.1,HEM.Bcin12g06390.1.col.0,HEM.Bcin12g06400.1.coi.1,HEM.Bcin12g06400.1.col.0,HEM.Bcin12g06410.1.coi.1,HEM.Bcin12g06410.1.col.0,HEM.Bcin12g06420.1.coi.1,HEM.Bcin12g06430.1.coi.1,HEM.Bcin12g06430.1.col.0)

#total sig SNPs per trait here
table(HEM.topSNPs$Trait)
jpeg(paste("plots/BOT_trueMAF20_NA10_lowTR_MAF20.ManhattanPlot.jpg", sep=""), width=7.5, height=5, units='in', res=600)
plot(ggplot(HEM.topSNPs, aes(x=Index, y=Effect))+
       theme_bw()+
       geom_point(aes(color = factor(Trait))))+
       geom_hline()
dev.off()
#greyscale version
#4 to 15
for (i in c(4:24)){
  #jpeg(paste("paper/plots/FigR6/bw_Sl_LesionSize_trueMAF20_NA10_lowTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=7.5, height=5, units='in', res=600)
  plot(ggplot(HEM.plotdata, aes(x=Index, y=100*HEM.plotdata[,i]))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Chrom)))+
         labs(list(y=expression(paste("Estimated Effect Size (",mm^{2},")")), title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
         theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         geom_hline(yintercept=100*get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
         geom_hline(yintercept=100*get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
         geom_text(aes(0,100*get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", size=14, vjust = 2, hjust=.05), col = "black")+
         geom_hline(yintercept=100*get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_hline(yintercept=100*get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_text(aes(0,100*get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", size=14, vjust = 4.5, hjust=.05), col = "black")+
         theme(legend.position="none")+
         theme(panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  )
  #dev.off()
}

#NA10 chromosomes
#scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+