#Nicole E Soltis
#

#plot of trans hotspots from GEMMA

#-------------------------------------------------------
rm(list=ls())
#setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_Bc")
mytop10 <- read.table("05_GEMMAsumm/GEMMA_top10SNPsample.txt", sep=",")
mytop100 <- read.table("05_GEMMAsumm/GEMMA_top100SNPsample.txt", sep=",")
mytop1 <- read.table("05_GEMMAsumm/GEMMA_top1SNPsample.txt", sep=",")
## check which SNP set here
mytop10 <- mytop1
mytop10$binfx <- 1

#read in full GEMMA output from 1 gene: set binary SNP fx to = 0. Assures plotting over whole genome, not just hotspot subsamples
allmySNP <- read.table("04_GEMMAoutput/binMAF20NA10_PLINK_1.assoc.txt")
names(allmySNP) <- c("chr", "rs", "ps", "n_mis", "n_obs","allele1", "allele0", "af", "beta" , "se", "p_wald", "p_lrt", "p_score")
allmySNP$pheno <- NA 
allmySNP$binfx <- 0

mydat <- rbind(allmySNP, mytop10)
which(grepl("chr", mydat$chr))
mydat <- mydat[-c(1),]
myGEMMA <- mydat

myGEMMA$chr <- as.numeric(myGEMMA$chr)
myGEMMA$ps <- as.numeric(myGEMMA$ps)
myGEMMA <- myGEMMA[with(myGEMMA, order(chr, ps)), ]

library("ggplot2")
#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20")
names(myColors) <- levels(myGEMMA$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#Make plotting variables
myGEMMA$Index = NA
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(myGEMMA$chr)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    myGEMMA[myGEMMA$chr==i, ]$Index=myGEMMA[myGEMMA$chr==i, ]$ps
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(myGEMMA,myGEMMA$chr==i-1)$ps, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    myGEMMA[myGEMMA$chr==i, ]$Index=myGEMMA[myGEMMA$chr==i, ]$ps+lastbase
  }
}

##check which df
#write.csv(myGEMMA, "05_GEMMAsumm/AllBcgenes_top100SNP_MAF20NA10_GEMMA_kmat1_Indexed.csv")

hist(myGEMMA$ps)
hist(myGEMMA$Index)
#positions look fine...

#summary plot: SNP overlap
mydat <- myGEMMA[,c("chr","ps","binfx","Index","pheno")]

library("dplyr")
mydat <- myGEMMA %>%
  select(chr, ps, binfx, Index) %>%
  group_by(Index) %>%
  #chr and ps should be identical per Index
  summarise(chr = mean(chr), ps = mean(ps), numsig = sum(binfx))

#num sig plot - currently by Index but could do by window
##check labeling for SNP number
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
jpeg(paste("plots/AllBcgenes_MAF20_10NA_GEMMA_top1SNP_hotspots.jpg", sep=""), width=8, height=5, units='in', res=600)
ggplot(mydat, aes(x=Index, y=numsig))+
  theme_bw()+
  colScale+
  geom_point(aes(color = factor(chr),alpha=0.001))+
  labs(list(y="Number of genes with a top 1 SNP", title="meta-plot across all 9k Bc genes"))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  theme(legend.position="none")+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #same for all 3 phenos
  scale_x_continuous(name="Chromosome", breaks = c(114174,  342554,  570966,  799392, 1028058, 1256043, 1484483, 1712912, 1941358, 2169675, 2398096, 2626508, 2854952, 3083366, 3311802, 3540153), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0)
dev.off()


setwd("~/Projects/BcAt_RNAGWAS")
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")
#jpeg(paste("plots/AtBclesion_MAF20_10NA_GEMMA_pad3.jpg", sep=""), width=8, height=5, units='in', res=600)
#print(ggplot(myGEMMA, aes(x=Index, y=beta))+
#columns with pscore to plot: 13 Col0, 16 coi1, 19 npr1, 22 pad3
#rows is pheno 1:4

print(
  ggplot(myGEMMA, aes(x=Index, y=binfx))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(chr),alpha=0.001))+
        labs(list(y="Top 10 SNP", title="meta-plot across all 9k Bc genes"))+
        guides(col = guide_legend(nrow = 8, title="Chromosome"))+
        theme(legend.position="none")+
        theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        #same for all 3 phenos
        scale_x_continuous(name="Chromosome", breaks = c(114174,  342554,  570966,  799392, 1028058, 1256043, 1484483, 1712912, 1941358, 2169675, 2398096, 2626508, 2854952, 3083366, 3311802, 3540153), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
        expand_limits(y=0)
)
dev.off()


#get chromosome midpoints
my.chroms <- as.data.frame(myGEMMA[!duplicated(myGEMMA$chr, fromLast=FALSE), "Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- myGEMMA[!duplicated(myGEMMA$chr, fromLast=TRUE), "Index"] # Upper Bounds
my.chroms$Chr.Mid <- round((my.chroms$Chr.Start + my.chroms$Chr.End)/2)
