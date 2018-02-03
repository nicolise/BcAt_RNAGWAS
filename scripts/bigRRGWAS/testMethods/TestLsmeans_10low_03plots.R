#reformat bigRR output data
#Nicole E Soltis
#013118

#--------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
setwd("~/Documents/GitRepos/BcAt_RNAGWAS")

library(tidyr); library(ggplot2); library(grid)

my.output <- read.csv("data/allreadsGWAS/03_bigRRout/testlowreads/lowreads_MAF20_013118.csv")
my.output <- my.output[,-c(1)]
my.output <- separate (my.output, outpt.HEM, into = c("Chrom", "Contig", "Pos") ) #Too many values at 3 locations: 5314, 79662, 217455
my.output$Chr.Con.Pos <- paste(my.output$Chrom, my.output$Contig, my.output$Pos, sep=".")
my.output$Chr.Con <- paste(my.output$Chrom, my.output$Contig, sep=".")

#sort dataframe rows in order of Chrom, then Pos
my.output <- my.output[with(my.output, order(Chrom, Contig, Pos)), ]
my.output <- my.output[,c(1:3,40,41,4:39)]
#now z-scale the phenotypes
my.output.z <- my.output[,c(1:5)]
for (i in c(6:41)){
  my.output.z[,i] <- scale(my.output[,i], center = TRUE, scale = TRUE)
  names(my.output.z)[i] <- paste(names(my.output)[i], "zeff", sep=".")
}
my.output.z <- my.output.z[,c(4,6:41)]
my.data <- merge(my.output,my.output.z, by="Chr.Con.Pos")

#----------------------------------------------------------------
#indexing:
#make sure things are ordered correctly
my.data$Chrom <- as.numeric(my.data$Chrom)
my.data$Contig <- as.numeric(my.data$Contig)
my.data$Pos <- as.numeric(my.data$Pos)
my.data <- my.data[with(my.data, order(Chrom, Contig, Pos)), ]
#let's try making the chrom.Cont integers so that R isn't confused
my.data$Chr.Con.F <- as.factor(my.data$Chr.Con)
unique(my.data$Chr.Con.F)
recode.vars <- data.frame(OGvals = as.factor(c("1.0", 1.1, "2.0", 2.1, 2.2, 2.3, 2.4, 2.5, "3.0", 3.1, 3.2, "4.0",  4.1, "5.0",  5.1, "6.0", 6.1, 6.2, 6.3, "7.0",  7.1, 7.2, "8.0",  8.1, 8.2, "9.0",  9.1, "10.0", 10.1, "11.0", 11.1, "12.0", 12.1, "13.0", 13.1, 13.2, "14.0", 14.1, 14.2, "15.0", 15.1, 15.2, 15.3, 15.4, "16.0", 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, "16.10", 16.11)), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56))

my.data$Chr.Con.Int <- recode.vars$newvals[match(my.data$Chr.Con.F, recode.vars$OGvals)]
unique(my.data$Chr.Con.Int)

#Make plotting variables -- a continuous count from Chromosome 1, Contig 1, Position 1 to the end of the last Contig of Chromosome 16.
my.data$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(my.data$Chr.Con.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    my.data[my.data$Chr.Con.Int==i, ]$Index=my.data[my.data$Chr.Con.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(my.data,my.data$Chr.Con.Int==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    my.data[my.data$Chr.Con.Int==i, ]$Index=my.data[my.data$Chr.Con.Int==i, ]$Pos+lastbase
  }
}

#--------------------------------------------------------------------
#plotting
#something is wrong in reshaped data frame: does not look anything like hand drawn plots
#first, scatterplots to compare effect sizes
name.key <- as.data.frame(names(my.data))

attach(my.data)
plot1 <- ggplot(data=my.data, aes(x=my.data[,6], y=my.data[,(6+12)]))+ geom_point(alpha=1/10)
plot2 <- ggplot(data=my.data, aes(x=my.data[,6], y=my.data[,(6+24)]))+ geom_point(alpha=1/10)

#plot for each transcript out of 12

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

#6:17
for (i in 6:17){
  #my.plot.dat <- splitlist[[i]]
  my.transcript <- names(my.data)[i]
  plot1 <- ggplot(data=my.data, aes(x=my.data[,i], y=my.data[,(i+12)]))+ geom_point(alpha=1/10) + labs(x=names(my.data)[i], y=names(my.data)[i+12]) + theme_bw()
  plot2 <- ggplot(data=my.data, aes(x=my.data[,i], y=my.data[,(i+24)]))+ geom_point(alpha=1/10) + labs(x=names(my.data)[i], y=names(my.data)[i+24]) + theme_bw()
  plot3 <- ggplot(data=my.data, aes(x=my.data[,(i+12)], y=my.data[,(i+24)]))+ geom_point(alpha=1/10) + labs(x=names(my.data)[i+12], y=names(my.data)[i+24]) + theme_bw()

# 2 figures arranged in 2 rows and 1 columns
file.name <- paste("plots/testMethods/scatter",my.transcript,".jpg", sep="")
jpeg(file.name, width=7.5, height=5, units='in', res=600)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(2, 1))
dev.off()
}

#-----------------------------------------------------------------------
#then, manhattan plots
# for finding the chromosome starts and finish for Chromosome Label Location:
my.chroms <- as.data.frame(my.data[!duplicated(my.data$Chrom, fromLast=FALSE), "Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- my.data[!duplicated(my.data$Chrom, fromLast=TRUE), "Index"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
names(my.data)
#6:77
for (i in c(6:77)){
  jpeg(paste("plots/testMethods/Manhattans/MAF20_", names(my.data)[i], ".ManhattanPlot.jpg", sep=""), width=7.5, height=5, units='in', res=600)
  plot(ggplot(my.data, aes(x=Index, y=my.data[,i]))+
         theme_bw()+
         #colScale+
         geom_point(aes(color = factor(Chrom)), alpha=1/2)+
         labs(list(y=expression(paste("Estimated Effect Size")), title=paste("Expression of ", names(my.data)[i])))+
         theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         theme(legend.position="none")+
         theme(panel.border = element_blank(), 
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))+
         scale_x_continuous(name="Chromosome", breaks = c(1677874,  5250045,  9006852, 11066684, 13584732, 17193234, 20022484, 22388798, 24412408, 26786090, 28588796, 30134140, 31893641, 34010154, 35809318, 38946222), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
         expand_limits(y=0))
  dev.off()
}



#-----------------------------------------------------------------------
#broken data reshaping here

#going to remove transcript ID info (just gene ID info)
colnames(my.data) <- gsub("\\.2\\.", "\\.1\\.", colnames(my.data))
colnames(my.data) <- gsub("\\.1\\.", "\\_", colnames(my.data))
attach(my.data)

births.long1<-reshape(births.wide, varying=c(“b2_01″,”b2_02″,”b2_03”, “b4_01”, “b4_02”, “b4_03″), direction=”long”, idvar=”caseid”, sep=”_”)

#go wide to long
my.data.long<-reshape(my.data, varying=c(6:77), direction="long", v.names=c("Bcin09g00820.1","Bcin11g00300.1","Bcin01g03900.1","Bcin01g03910.1","Bcin01g03910.2","Bcin16g05200.1","Bcin09g01130.1","Bcin10g00140.1","Bcin10g00010.1","Bcin14g05180.1","Bcin07g04990.1","Bcin10g00150.1"), times=c("HEM","cutoff.HEM","zscale.HEM","HEM.zeff","cutoff.HEM.zeff","zscale.HEM.zeff"))

my.data.long.b<-reshape(my.data, varying=c(6:77), direction="long", times=c("Bcin09g00820.1","Bcin11g00300.1","Bcin01g03900.1","Bcin01g03910.1","Bcin01g03910.2","Bcin16g05200.1","Bcin09g01130.1","Bcin10g00140.1","Bcin10g00010.1","Bcin14g05180.1","Bcin07g04990.1","Bcin10g00150.1"), v.names=c("HEM","cutoff.HEM","zscale.HEM","HEM.zeff","cutoff.HEM.zeff","zscale.HEM.zeff"))

names(my.data.long)[9]<- "Method"
names(my.data.long.b)[9] <- "Transcript"

#first, split by transcript

#then, plots for each!
splitlist <- split(my.data.long.b, my.data.long.b$Transcript)

for (i in 1:12){
  my.plot.dat <- splitlist[[i]]
  my.transcript <- rownames(splitlist[[i]])[1]
  plot1 <- ggplot(data=my.plot.dat, aes(x=HEM, y=cutoff.HEM))+ geom_point(alpha=1/10)
  plot2 <- ggplot(data=my.plot.dat, aes(x=HEM, y=zscale.HEM))+ geom_point(alpha=1/10)
  plot3 <- ggplot(data=my.plot.dat, aes(x=cutoff.HEM, y=zscale.HEM))+ geom_point(alpha=1/10)
  
  # 2 figures arranged in 2 rows and 1 columns
  file.name <- paste("plots/scatter",my.transcript,".jpg", sep="")
  jpeg(file.name, width=7.5, height=5, units='in', res=600)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(plot1, vp = vplayout(1, 1))
  print(plot2, vp = vplayout(1, 2))
  print(plot3, vp = vplayout(2, 1))
  dev.off()
}