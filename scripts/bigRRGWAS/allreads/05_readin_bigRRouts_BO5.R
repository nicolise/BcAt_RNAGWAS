#Nicole E Soltis
#020518
#---------------------------------------------------------
library(tidyr);library(ggplot2)
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/data/allreadsGWAS/BO5.10/03_bigRRout/subset/")
#read in individual bigRR output files (1 per geno)
#read in all files in folder by pattern matching

#testing: just files on Chr1

my.files <- list.files(pattern = ".csv")
pList <- list()
for(i in 1:length(my.files)) {
  my.file <- read.csv(my.files[i])
  my.file <- my.file[,-c(1)]
  pList[[i]] <- my.file
}

#merge all phenotypes into one file -- may only make sense for small sets
my.full.data <- pList[[1]]
for(j in 2:length(pList)){
  my.full.data <- cbind(my.full.data, pList[[j]][,2])
  names(my.full.data)[j+1] <- names(pList[[j]])[2]
}

#so: do top SNP locations correlate with gene location? (Cis fx)
#need: key for gene locations from Vivian
#first step: plot top 100 SNPs per transcript
my.full.data<-my.full.data[,order(colnames(my.full.data),decreasing=TRUE)]
names(my.full.data)

#get gff3 file info
library("ape")
setwd("~/Projects/BcGenome/data/ensembl/BO5.10")
#can extract the files from .gz using 7-zip
#then read using R!
my.gff <- read.gff("extracted/Botrytis_cinerea.ASM83294v1.38.chromosome.1.gff3/Botrytis_cinerea.ASM83294v1.38.chromosome.1.gff3", na.strings = c(".", "?"))
my.gene.list <- as.data.frame(names(my.full.data)[2:31])
names(my.gene.list)[1] <- "mygenes"
my.gene.list$mygenes <- as.character(my.gene.list$mygenes)
#some lazy regex
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$mygenes, ".1.HEM"), "[", 1)
my.gene.list$justgenes <- sapply(strsplit(my.gene.list$justgenes, ".2.HEM"), "[", 1)
#now loop over to keep things
#1:length(my.gene.list$justgenes)
my.gff.genes <- my.gff[1,]
my.gff.genes$transcript <- NA
my.gff.genes <- my.gff.genes[-c(1),]
for (y in (1:length(my.gene.list$justgenes))){
  my.current.gene <- my.gff[grep(my.gene.list[y,2], my.gff$attributes), ]
  my.current.gene$transcript <- my.gene.list[y,2]
  my.gff.genes <- rbind(my.gff.genes, my.current.gene)
}

#for a first go, just keep the "gene" rows
my.transcript.locs <- my.gff.genes[my.gff.genes$type=="gene", ]

#---------------------------------------------------------------------------------
#now: plot for each transcript
#my.full.data includes effect sizes for each transcript across whole genome
#first plot: manhattan of each transcript individually, with a bar for gene location
my.full.data <- separate(my.full.data, outpt.HEM, into = c("Chrom","Pos") ) 
my.full.data$Chr.Pos <- paste(my.full.data$Chrom, my.full.data$Pos, sep=".")

#sort dataframe rows in order of Chrom, then Pos
my.plotdat <- my.full.data
my.plotdat$Chrom <- as.numeric(my.plotdat$Chrom)
my.plotdat$Pos <- as.numeric(my.plotdat$Pos)
my.plotdat<- my.plotdat[with(my.plotdat, order(Chrom, Pos)), ]

#Add indexing
#Make plotting variables -- a continuous count from Chromosome 1, Contig 1, Position 1 to the end of the last Contig of Chromosome 16.
my.plotdat$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(my.plotdat$Chrom)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    my.plotdat[my.plotdat$Chrom==i, ]$Index=my.plotdat[my.plotdat$Chrom==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(my.plotdat,my.plotdat$Chrom==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    my.plotdat[my.plotdat$Chrom==i, ]$Index=my.plotdat[my.plotdat$Chrom==i, ]$Pos+lastbase
  }
}

#then, manhattan plots
# for finding the chromosome starts and finish for Chromosome Label Location:
my.chroms <- as.data.frame(my.plotdat[!duplicated(my.plotdat$Chrom, fromLast=FALSE), "Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- my.plotdat[!duplicated(my.plotdat$Chrom, fromLast=TRUE), "Index"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2

#3:32
#whole genome
for (j in c(3:32)){
  jpeg(paste("plots/BO5/Manhattans/cis_test/MAF20_", names(my.plotdat)[j], ".ManhattanPlot.jpg", sep=""), width=7.5, height=5, units='in', res=600)
  plot(ggplot(my.plotdat, aes(x=Index, y=my.plotdat[,j]))+
         theme_bw()+
         #colScale+
         geom_point(aes(color = factor(Chrom)), alpha=1/2)+
         labs(list(y=expression(paste("Estimated Effect Size")), title=paste("Expression of ", names(my.plotdat)[j])))+
         theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         theme(legend.position="none")+
         theme(panel.border = element_blank(), 
               #panel.grid.major = element_blank(),
               #panel.grid.minor = element_blank(), 
               axis.line = element_line(colour = "black"))+
         scale_x_continuous(name="Chromosome", breaks = c(2043070, 5756531, 9040102, 11884071, 14583781, 17407831, 20088809, 22716031, 25287492, 27760698, 30133998, 32479787, 34784035, 36981647, 39063874, 41056728, 42151808, 42364627), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))+
         expand_limits(y=0))
  dev.off()
}

#just chromosome 1, with gene location plotted
my.plotchr1 <- my.plotdat[my.plotdat$Chrom==1,]
#3:32
for (j in c(3:32)){
  my.get.gene <- names(my.plotchr1)[j]
  #then look that name up in my.gene.list under mygenes
  my.get.transcript <- my.gene.list[my.gene.list$mygenes==my.get.gene,2]
  #then look up the corresponding justgenes in my.transcript.locs
  #then extract gene start and end
  my.gene.start <- my.transcript.locs[my.transcript.locs$transcript==my.get.transcript,4]
  my.gene.end <- my.transcript.locs[my.transcript.locs$transcript==my.get.transcript,5]
  jpeg(paste("plots/BO5/Manhattans/cis_test/MAF20_C1_", names(my.plotchr1)[j], ".ManhattanPlot.jpg", sep=""), width=7.5, height=5, units='in', res=600)
  plot(ggplot(my.plotchr1, aes(x=Index, y=my.plotchr1[,j]))+
         theme_bw()+
         #colScale+
         geom_point(aes(color = factor(Chrom)), alpha=1/2)+
         labs(list(y=expression(paste("Estimated Effect Size")), title=paste("Expression of ", names(my.plotchr1)[j])))+
         theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         theme(legend.position="none")+
         theme(panel.border = element_blank(), 
               #panel.grid.major = element_blank(),
               #panel.grid.minor = element_blank(), 
               axis.line = element_line(colour = "black"))+
         #gene start stop 883036 - 883939 (positions on chromosome 1) plus 10kb each side
         geom_rect(mapping=aes(ymin=(0.8*max(my.plotchr1[,j])), ymax=(max(my.plotchr1[,j])), xmin=(my.gene.start-10000), xmax=(my.gene.end+10000)), alpha=0.01, fill="darkturquoise")+
         expand_limits(y=0))
  dev.off()
}
