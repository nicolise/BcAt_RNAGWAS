#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())

#rename files
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_les/troubleshoot_Col0/D_04_GEMMAout/")
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_les/troubleshoot_Col0/D_04_GEMMAout/")
#read in each of 12 phenotypes
my.files <- list.files(pattern = c("assoc"))
# #2_LA1589 fails with both k-matrix options.
# my.names <- c("1_Col0.Les","2_coi1.Les","3_npr1.Les","4_pad3.Les","5_tga3.Les","6_anac055.Les")

# #rename all files
# for(i in 1:length(my.files)) {
#   my.file <- read.csv(my.files[i])
#   file.rename(from=file.path(my.files[i]), to=file.path(paste(file.path(my.files[i]),my.names[i],".txt",sep="")))
# }

# my.files <- list.files(pattern = c("assoc"))
# for (i in 1:length(my.files)){
#   my.file <- read.table(my.files[i], header=TRUE)
#   my.file$chr.ps <- paste(my.file$chr, my.file$ps, sep=".")
#   ifelse(i == 1, full.file <- my.file, full.file <- merge(full.file, my.file[,c("beta","se","p_score","chr.ps")], by="chr.ps"))
#   ifelse(i == 1, names(full.file)[9] <- paste(my.names[i], "_beta", sep=""), names(full.file)[(ncol(full.file)-2)] <- paste(my.names[i], "_beta", sep=""))
#   ifelse(i == 1, names(full.file)[10] <- paste(my.names[i], "_se", sep=""), names(full.file)[(ncol(full.file)-1)] <- paste(my.names[i], "_se", sep=""))
#   ifelse(i == 1, names(full.file)[13] <- paste(my.names[i], "_pscore", sep=""), names(full.file)[(ncol(full.file))] <- paste(my.names[i], "_pscore", sep=""))
# }
# full.file <- full.file[,-c(12,13)]

setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_les")
#write.csv(full.file, "D_06_results/LesionPhenos_allSNPs_MAF20NA10_GEMMA_kmat1.csv")
full.file <- read.csv("D_06_results/LesionPhenos_allSNPs_MAF20NA10_GEMMA_kmat1.csv")
library(ggplot2); 

#let's try a manhattan plot. Choosing score test for now.

#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(full.file$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#sort dataframe rows in order of Chrom, then Pos
str(full.file)

##test: read pheno 1 directly as full.file, see if it fixes plot
#full.file <- read.table("D_04_ogphenos/binMAF20NA10_norand_kmat1_pheno1.assoc.txt1_Col0.Les.txt", header=TRUE)

myGEMMA <- full.file[with(full.file, order(chr, ps)), ]

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

#comparing positions to B05.10 outputs in 06_transPatterns_B05.R
#pos should be about 0 to 14,000
#ps goes 121 to 408,6062 ... is one of the chromosomes too long somehow?
#Index should be about 1000 to 41,997,650
#Index goes 4224 to 40,421,038
hist(myGEMMA$ps)
hist(myGEMMA$Index)
#positions look fine...

#get thresholds here 
mythrs <- #read.csv("data/GEMMA_files/D_07_randOUTS/GEMMA_1krand_thresholds.csv")
mythrs

#troubleshooting col0
col0pval <- as.data.frame(table(myGEMMA$X1_Col0.Les_pscore))
View(col0pval)
topcol0pval <- col0pval[col0pval$Freq %in% c(408,917:6075),]
lowcol0pval <- col0pval[!col0pval$Freq %in% c(408,917:6075),]
editcol0 <- myGEMMA[!myGEMMA$X1_Col0.Les_pscore %in% topcol0pval$Var,]
# basically looks like col0 failed. minimum "believable" p-val is 0.0396

setwd("~/Projects/BcAt_RNAGWAS")
jpeg(paste("plots/AtBclesion_MAF20_10NA_GEMMArand_col0.jpg", sep=""), width=8, height=5, units='in', res=600)
#print(ggplot(myGEMMA, aes(x=Index, y=beta))+
#columns with pscore to plot: 12, 15, 18, 21
##print(ggplot(myGEMMA, aes(x=Index, y=(-log10(myGEMMA[,13]))))+
ggplot(editcol0, aes(x=Index, y=(-log10(editcol0$X1_Col0.Les_pscore))))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(chr),alpha=0.001))+
        labs(list(y=expression('-log'[10]*'p'), title="Col-0 Lesion Size"))+
        guides(col = guide_legend(nrow = 8, title="Chromosome"))+
        #geom_hline(yintercept=-log10(), colour = "black", lty=2)+ #250
        #geom_hline(yintercept=-log10(), colour = "black", lty=3)+ #2500
        theme(legend.position="none")+
        theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        #same for all 3 phenos
        scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
        expand_limits(y=0)
)
dev.off()


#get chromosome midpoints
my.chroms <- as.data.frame(myGEMMA[!duplicated(myGEMMA$chr, fromLast=FALSE), "Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- myGEMMA[!duplicated(myGEMMA$chr, fromLast=TRUE), "Index"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2

#now read in all 3, make combination file for meta-analysis