#Nicole E Soltis
#08/

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#At transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At/output
#kmat script: 
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At/A04_runGEMMA_kmat_3genos.sh

rm(list=ls())
setwd("/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At")

#read in individual GEMMA output files (1 per geno)
#read in all files in folder by pattern matching
#my.files <- list.files(pattern = ".assoc.txt")
    
  #each phenotype
for (j in c("col0")){
 for (i in 1:23947){
    #actually: 1:23947 for At
    Sys.time()
    my_gemma <- read.table(paste("output/",j,"_MAF20NA10_obs_",i,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 4 seconds to read 1 phenotype
    #times 24000 = 27 hours
    #take top 10 SNP/phenotype
    #also save top 1 SNP/ phenotype
    my_gemma$pheno <- i
    my_gemma_top100 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top100 <- my_gemma_top100[1:100,]
    row.names(my_gemma_top100) <- c(1:100) + (100*(i-1))
    my_gemma_top10 <- my_gemma_top100[1:10,]
    my_gemma_top1 <- my_gemma_top10[1,]
    #and z scaling
    my_gemma.z <- my_gemma
    my_gemma.z$beta_z <- scale(my_gemma.z$beta, center = TRUE, scale = TRUE)
    my_gemma.z <- my_gemma.z[abs(my_gemma.z$beta_z) > 4,]
    mylgsnp <- my_gemma[abs(my_gemma$beta) > 0.5,]
    #this gives an error but it's fine
    try(ifelse( i == 1, write.table(my_gemma_top100, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top100SNPsample.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma_top100, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top100SNPsample.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top10, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top10SNPsample.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma_top10, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top10SNPsample.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top1, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top1SNPsample.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma_top1, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top1SNPsample.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma.z, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_topSNPsample_zscale.txt", sep=""), sep = ",", col.names = TRUE), write.table(my_gemma.z, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_topSNPsample_zscale.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(mylgsnp, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top100_beta05SNP.txt", sep=""), sep = ",", col.names = TRUE), write.table(mylgsnp, paste("05_GEMMAsumm/",j,"_round2/",j,"_GEMMA_top100_beta05SNP.txt", sep=""), sep = ",", col.names = FALSE, append = TRUE)))
    Sys.time()
  }
}
#-----------------------------------------------------------------------------------