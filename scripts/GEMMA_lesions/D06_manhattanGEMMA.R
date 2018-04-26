#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
#first round: just one file at a time. Then convert to loops
#13 is Domest, 14 is Wild, 15 is Sensitivity

#on laptop:
myGEMMA <- read.table("data/B05_GEMMA_les/D_04_ogphenos/binMAF20NA10_norand_kmat1_pheno1.assoc.txt", header=TRUE)

#rename files
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/B05_GEMMA_les/D_04_ogphenos/")
setwd("~/Projects/BcAt_RNAGWAS/data/B05_GEMMA_les/D_04_ogphenos/")
#read in each of 12 phenotypes
my.files <- list.files(pattern = c("assoc"))
#2_LA1589 fails with both k-matrix options.
my.names <- c("1_Col0.Les","2_coi1.Les","3_npr1.Les","4_pad3.Les","5_tga3.Les","6_anac055.Les")
#  "LA1547","LA1589","LA1684","LA2093","LA2176","LA2706","LA3008","LA3475","LA410","LA4345","LA4355","LA480"
#rename all files
# for(i in 1:length(my.files)) {
#   my.file <- read.csv(my.files[i])
#   file.rename(from=file.path(my.files[i]), to=file.path(paste(file.path(my.files[i]),my.names[i],".txt",sep="")))
# }

library(ggplot2); 

#virtually identical, just checking.
plot(myGEMMA$p_wald, myGEMMA$p_score)
plot(myGEMMA$p_wald, myGEMMA$p_lrt)
#wald test, likelihood ratio test, or score test

#let's try a manhattan plot. Choosing score test for now.

#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(myGEMMA$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#sort dataframe rows in order of Chrom, then Pos
str(myGEMMA)
#myGEMMA$ps <- as.numeric(myGEMMA$ps)
#myGEMMA$chr <- as.numeric(myGEMMA$chr)
myGEMMA <- myGEMMA[with(myGEMMA, order(chr, ps)), ]

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
mythrs <- read.csv("data/GEMMA_files/D_07_randOUTS/GEMMA_1krand_thresholds.csv")
mythrs

jpeg(paste("paper/plots/addGEMMA/SlBc_MAF20_10NA_GEMMArand_Sens.jpg", sep=""), width=8, height=5, units='in', res=600)
#print(ggplot(myGEMMA, aes(x=Index, y=beta))+
print(ggplot(myGEMMA, aes(x=Index, y=(-log10(p_score))))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(chr),alpha=0.001))+
        labs(list(y=expression('-log'[10]*'p'), title="Sensitivity"))+
        guides(col = guide_legend(nrow = 8, title="Chromosome"))+
        geom_hline(yintercept=-log10(2.350535e-03), colour = "black", lty=2)+ #250
        geom_hline(yintercept=-log10(1.419165e-02), colour = "black", lty=3)+ #2500
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