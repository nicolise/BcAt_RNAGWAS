#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")
#on windows laptop
setwd("~/Projects/BcSolGWAS/data")
#convert .tab SNP file to .csv
#sticking to MAF20, 10% missingness

 tab20 = read.delim("GWAS_files/01_tabfiles/Suzi_033016/Haploid_SNPS_97_dp6_maf20.tab")
  write.table(tab20, file="GEMMA_files/02_csvPrep/snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)


#convert each file to binary
SNPsMAF20 <- read.csv("GEMMA_files/02_csvPrep/snps_maf20.csv")
mySNPs <- SNPsMAF20

mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="."]<-NA #this is a true NA
allSNPs<- mySNPs

#remove low MAFs!
names(mySNPs)

mySNPs$Freq.A <- rowSums(mySNPs[,c(4:100)] =="A", na.rm=T)
mySNPs$Freq.C <- rowSums(mySNPs[,c(4:100)] =="C", na.rm=T)
mySNPs$Freq.T <- rowSums(mySNPs[,c(4:100)] =="T", na.rm=T)
mySNPs$Freq.G <- rowSums(mySNPs[,c(4:100)] =="G", na.rm=T)
mySNPs$All <- (mySNPs$Freq.A + mySNPs$Freq.C + mySNPs$Freq.G + mySNPs$Freq.T)
#97 isolates total
#remove all loci (rows) with >10% missingness (mySNPs$All <87.3)
mySNPs <- mySNPs[mySNPs$All>87,] 

#now remove all loci (rows) with MAF <20 (max > 77.6)
mySNPs <- transform(mySNPs, max = pmax(Freq.A, Freq.C, Freq.G, Freq.T))
mySNPs <- mySNPs[mySNPs$max < 78,] #this is fine, keep all.
mySNPs$missing <- apply(mySNPs[,c(4:100)], 1, function(x) sum(is.na(x)))

#MAF is too low if: max + missing > 78
mySNPs <- mySNPs[(mySNPs$max + mySNPs$missing) < 78,]

#and now for making PED format for PLINK!
  #do not need positional info: just SNP states for PED
#turn df sideways (individuals as rows, SNPs as columns)
#split each genotype into 2 identical columns (PED assumes diploid)
#add a first column: FAM1 (no info on isolate families)
#second column: isolate ID
#third column: father ID (a column of zeros)
#fourth column: mother ID (a column of zeros)
#fifth column: individual sex = 1 (all assumed same)
#sixth  column: binary  phenotype (all = 1)
#fix column order
mySNPs2 <- mySNPs[,-c(1:3, 101:107)]
#write.csv(mySNPs2, "GEMMA_files/02_csvPrep/hp_charMAF20_10NA_forPED.csv")
#df2[,c(1,3,2,4)]

#in 48 + 19 hours, finished 465,853 rows = 85% of 546,540
#I'll finish from there. 465853/2 = 232926 is done. 232927 was half done
#mySNPs3 <- rbind(mySNPs3, mySNPs2[232927,])
#now start from 232928 
#second time through: [1] "2017-07-05 13:52:59 PDT" until [1] "2017-07-09 14:22:10 PDT", so 4 days

# Sys.time()
# mySNPs3 <- as.data.frame(NULL)
# for(i in 1:nrow(mySNPs2)){
#   mySNPs3<-rbind(mySNPs3, mySNPs2[i,])
#   mySNPs3<-rbind(mySNPs3, mySNPs2[i,])
# }
# Sys.time()
write.csv(mySNPs3, "GEMMA_files/02_csvPrep/dp_charMAF20_10NA.csv")

mySNPs3 <- read.csv("GEMMA_files/02_csvPrep/fulldata/dp_charMAF20_10NA.csv")
#transpose and format for PED
mySNPs3 <- mySNPs3[,-c(1)]
mySNPs4 <- as.data.frame(t(mySNPs3))
#add binary phenotype = 1 (6)
mySNPs4 <- cbind("Pheno" = 1, mySNPs4)
#add individual sex = 1 (5)
mySNPs4 <- cbind("sex" = 1, mySNPs4)
#add Mother = 0 (4)
mySNPs4 <- cbind("Mother" = 0, mySNPs4)
#add Father = 0 (3)
mySNPs4 <- cbind("Father" = 0, mySNPs4)
#turn row names into column 2
mySNPs4 <- cbind(rownames(mySNPs4), mySNPs4)
colnames(mySNPs4)[1] <- 'Isolate'
#add the fam column (1)
mySNPs4 <- cbind("FAM" = "FAM1", mySNPs4)

#make a MAP file for plink (need it to make the bed (binary ped) file from ped)
myMAP <- mySNPs[,c("X.CHROM","POS")]
#remove "Chromosome" from X.CHROM
myMAP2 <- as.data.frame(lapply(myMAP, function(x) {
                 gsub("Chromosome", "", x)
              }))
write.csv(myMAP2, "GEMMA_files/02_csvPrep/fulldata/MAP_dpcharMAF20NA10.csv")
#add a column of "SNP identifiers" in excel and remove headers

write.csv(mySNPs3, "GEMMA_files/02_csvPrep/fulldata/dp_charMAF20_10NA.csv")
write.csv(mySNPs, "GEMMA_files/02_csvPrep/fulldata/hp_charMAF20_10NA.csv")
write.csv(mySNPs4, "GEMMA_files/02_csvPrep/fulldata/PED_dpcharMAF20NA10.csv")
