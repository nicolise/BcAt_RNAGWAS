#Nicole E Soltis
#04/05/18

#-------------------------------------------------------------------------
rm(list = ls())
#list all files in /media/ outputs
#loop through all 1000 files
# for each of 12 phenotypes, create new dataframe that keeps only selected SNPs:
  #df 1: ordered on SNP p-value top 1% = 2720 SNPs ish
  #df 2: ordered on SNP p-value small -> large, take 1st, 2nd, 25th, 250th, 2500th SNP
#save df 1 to media
#save df 2 to ~/Projects
#loop through all df2 to build a new df across all 1000 permutations using rbind
  #columns: chromosome, position, p-value on phenotype 1...12, which SNP (1 / 2 / 25 / 250 / 2500)
  #
#setwd("~/Projects/BcSolGWAS/data/GEMMA_files/D_05_bigrand")
setwd("/media/nesoltis/Data/Kliebenstein/Soltis/BcSolGWAS/data/GEMMA_files/")

#read in files from each of 1000 rand folders
#each randomization
#for (i in 1:1000){
#done: i = 1:500
#60 phenos in 2 hours
#started 561:1000 at 3:20pm 04/10/18 ... expect 15 hours to complete = 6:30am
for (i in 561:1000){
  newdir <- paste0("D_06_randOUT/topouts/rand1k_",i)
  dir.create(newdir)
  #each phenotype
  for (j in 1:15){
    Sys.time()
    my_gemma <- read.table(paste("D_06_randOUT/output/rand1k_",i,"/pheno",j,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 8 seconds to read 1 phenotype
    #times 15 times 1000 = 33 hours
    ntop1 <- nrow(my_gemma)*0.01
    my_gemma_top1 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top1 <- my_gemma_top1[1:2500,]
    my_gemma_points <- my_gemma_top1[c(1,2,25,250,2500),]
    my_gemma_points$SNPnum <- c(1,2,25,250,2500)
    my_gemma_points$run <- i
    my_gemma_points$pheno <- j
    #this gives an error but it's fine
    try(ifelse( i == 1 & j ==1, write.table(my_gemma_points, "D_07_randSUMM/GEMMA_1krand_SNPsample.csv", sep = ",", col.names = T), write.table(my_gemma_points, "D_07_randSUMM/GEMMA_1krand_SNPsample.csv", sep = ",", col.names = F, append = T)))
    write.csv(my_gemma_top1, paste("D_06_randOUT/topouts/rand1k_",i,"/pheno",j,".csv", sep=""))
    Sys.time()
  }
}

#now cp D_07_randSUMM to ~/Documents/GitRepos/BcSolGWAS/
#and make sure nesoltis user account has rwx permission for the new files

#----------------------------------------------------------------
#now extract percentiles for QQ plots in each phenotype!
setwd("/media/nesoltis/Data/Kliebenstein/Soltis/BcSolGWAS/data/GEMMA_files/")

for (i in 1:1000){
  #each phenotype
newdir <- paste0("D_06_randOUT/quantiles/rand1k_",i)
dir.create(newdir)
  for (j in 1:15){
    Sys.time()
    my_gemma <- read.table(paste("D_06_randOUT/output/rand1k_",i,"/pheno",j,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 8 seconds to read 1 phenotype
    #times 15 times 1000 = 33 hours
    
    #need to extract
    myquantsdf <- as.data.frame(NULL)
    for (myquant in seq(0.01,1,0.01)){
      quantout <- quantile(my_gemma$p_score, myquant)
      #print(quantout)
      myrow <- myquant*100
      myquantsdf[myrow,1] <- myrow
      names(myquantsdf)[1] <- "Quantile"
      myquantsdf[myrow,2] <- quantout
      names(myquantsdf)[2] <- "p_score"
      myquantsdf[myrow,3] <- i
      names(myquantsdf)[3] <- "permutrun"
      myquantsdf[myrow,4] <- j
      names(myquantsdf)[4] <- "pheno"
    }
    try(ifelse( i == 1 & j ==1, write.table(myquantsdf, "D_07_randSUMM/GEMMA_1krand_SNPquants.csv", sep = ",", col.names = T), write.table(myquantsdf, "D_07_randSUMM/GEMMA_1krand_SNPquants.csv", sep = ",", col.names = F, append = T)))
     write.csv(myquantsdf, paste("D_06_randOUT/quantiles/rand1k_",i,"/pheno",j,".csv", sep=""))
    Sys.time()
  }
}
