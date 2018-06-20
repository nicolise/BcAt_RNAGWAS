#Nicole E Soltis
#06/08/18
#03_bigRR_cistrans
#--------------------------------------------------------------------
#extract BotBoaNet5 bigRR data to subfolder
# try post-hoc meta-analysis across phenotypes
  #first approach: MANTEL in perl? beta-value scaling
#check cis vs. trans SNP effect estimates
#later repeat this for GEMMA
rm(list = ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/")
#here are the original phenotypes for these SNPs
#total of 30 genes
PhenosNet <- read.csv("data/BcBotGWAS/02_MatchGenos/BOTBOANet5phenos.csv")

setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreads_bigRR/B05.10")
getPhenos <- as.data.frame(names(PhenosNet))
names(getPhenos)[1] <- "gene"
getPhenos$FileName <- paste("03_bigRRout_partial/outfiles/",getPhenos$gene, ".HEM.csv", sep="")
getPhenos <- getPhenos[-c(1:2),]

#file.copy(from=getPhenos$FileName, to="04_NetSubset/", overwrite = TRUE, recursive = FALSE, copy.mode = TRUE) #got 19/ 30. Now unzipping full files to find last 11

getPhenos$ZipFile <- paste("03_bigRRout/03_bigRRout/outfiles/",getPhenos$gene, ".HEM.csv", sep="")
#file.copy(from=getPhenos$ZipFile, to="04_NetSubset/", 
#          overwrite = TRUE, recursive = FALSE, 
#          copy.mode = TRUE) #all but one copied
getPhenos$NewFiles <- paste(getPhenos$gene, ".HEM.csv", sep="")

#now extract relevant GWAS data from here...
#from data/allreads_bigRR/B05.10/04_NetSubset/
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data/allreads_bigRR/B05.10/04_NetSubset")
my.files <- list.files(pattern = c(".HEM.csv"))
#for (i in 1:length(my.files)){
i <- 1
   my.file <- read.csv(my.files[i], header=TRUE)
   my.file$chr.ps <- paste(my.file$chr, my.file$ps, sep=".")
   ifelse(i == 1, full.file <- my.file, full.file <- merge(full.file, my.file[,c("beta","se","p_score","chr.ps")], by="chr.ps"))
   ifelse(i == 1, names(full.file)[9] <- paste(my.names[i], "_beta", sep=""), names(full.file)[(ncol(full.file)-2)] <- paste(my.names[i], "_beta", sep=""))
   ifelse(i == 1, names(full.file)[10] <- paste(my.names[i], "_se", sep=""), names(full.file)[(ncol(full.file)-1)] <- paste(my.names[i], "_se", sep=""))
   ifelse(i == 1, names(full.file)[13] <- paste(my.names[i], "_pscore", sep=""), names(full.file)[(ncol(full.file))] <- paste(my.names[i], "_pscore", sep=""))
 }
 full.file <- full.file[,-c(12,13)]

#read and combine all files


#next: check for haplotype structure?? haplotype-based model to locate cis effects