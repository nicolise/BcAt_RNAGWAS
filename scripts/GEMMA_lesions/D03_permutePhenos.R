#Nicole E Soltis
#R loop to run to generate randomized phenotypes, then feed into GEMMA
#-------------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")

#pipeline: see D02_prepPhenos.R
#have already run GEMMA for unpermuted phenotypes
#here: looped permutations to generate randomized phenotype files
#then: looped GEMMA across all files
#then: looped value extractions

#now randomize each phenotype! Wheee
Phenos_match <- read.csv("B05_GEMMA_les/D_02_randGEMMA/binMAF20NA10_fam.csv")

#to do this, open Rstudio using sudo in command line!
#need to save these on /media/ because there is not enough space on C://
#test on a small run 
#make directories on /media/ : BcSolGWAS/data/GEMMA_files/D_05_bigrand/
#and GEMMA_files/D_02_randGEMMA
#paste .bed and .bim into D_02_randGEMMA
setwd("/media/nesoltis/Data/Kliebenstein/Soltis/BcSolGWAS/data")
Sys.time()
for (i in 1:1000){
  Phenos_rand <- transform(Phenos_match, LA0410 = sample(LA0410), LA0480 = sample(LA0480), LA1547 = sample(LA1547), LA1589 = sample(LA1589), LA1684 = sample(LA1684), LA2093 = sample(LA2093), LA2176 = sample(LA2176), LA2706 = sample(LA2706), LA3008 = sample(LA3008), LA3475 = sample(LA3475), LA4345 = sample(LA4345), LA4355 = sample(LA4355),  Domesticated = sample(Domesticated), Wild = sample(Wild), DmWoD = sample(DmWoD))
  #select columns
  Phenos_rand <- Phenos_rand[,-c(1)]
  newdir <- paste0("GEMMA_files/D_05_bigrand/rand1k_",i)
  dir.create(newdir)
  cwd <- getwd()
  setwd(newdir)
  write.table(Phenos_rand, "binMAF20NA10_rand.fam", row.names=FALSE, col.names=TRUE)
  setwd(cwd)
  #and now copy .bed and .bim over to new directory
  plink.folder <- paste0(cwd,"/GEMMA_files/D_02_randGEMMA")
  new.folder <- paste0(cwd,"/",newdir)
  # find the files that you want
  my.bed <- list.files(plink.folder, ".bed$",full.names=T)
  my.bim <- list.files(plink.folder, ".bim$",full.names=T)
  # copy the files to the new folder
  file.copy(my.bed, new.folder)
  file.copy(my.bim, new.folder)
}
Sys.time()