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
#somehow 217749 of the non-duplicated original SNPs have become 217749 SNPs with 205 duplicated. Estimates for duplicated SNPs are very different --> going to drop all of these

full.file <- NULL
for (i in 1:length(my.files)){
   my.file <- read.csv(my.files[i], header=TRUE)
   names(my.file)[2] <- "chr.ps"
   print(sum(duplicated(my.file[,2])))
   my.file <- my.file[!duplicated(my.file$chr.ps),]
   ifelse(i == 1, full.file <- my.file, full.file <- merge(full.file, my.file[,c(2,3)], by="chr.ps"))
   #ifelse(i == 1, names(full.file)[9] <- paste(my.names[i], "_beta", sep=""), names(full.file)[(ncol(full.file)-2)] <- paste(my.names[i], "_beta", sep=""))
}
#10 mins to run
Sys.time()
write.csv(full.file, "BotBoaNet_allGenes_beta.csv")

#next: check for haplotype structure?? haplotype-based model to locate cis effects