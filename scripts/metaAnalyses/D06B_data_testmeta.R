#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#get thresholds here 
mythrs <- read.csv("data/GEMMA_files/D_07_randOUTS/GEMMA_1krand_thresholds.csv")

#LA2093 is pheno 6 (see Phenos_match)
Phenos_match <- read.csv("data/GEMMA_files/D_02_randGEMMA/binMAF20NA10_fam.csv")
names(Phenos_match)

#now read in all 12 phenos
setwd("~/Projects/BcSolGWAS/data/GEMMA_files/D_04_randphenos/")
        
  #read in each of 12 phenotypes
  # my.files <- list.files(pattern = c("assoc"))
  
  #check these: 
  names(Phenos_match)
  my.names <- c("1_LA0410", "10_LA3475", "11_LA4345", "12_LA4355", "13_Domest", "14_Wild", "15_Sens", "2_LA0480", "3_LA1547", "4_LA1589","5_LA1684","6_LA2093","7_LA2176","8_LA2706","9_LA3008")
  #rename all files
  # for(i in 1:length(my.files)) {
  #   my.file <- read.csv(my.files[i])
  #   file.rename(from=file.path(my.files[i]), to=file.path(paste(file.path(my.files[i]),my.names[i],".txt",sep="")))
  # }
  
  #now read in files and combine
  my.files <- list.files(pattern = c("assoc"))
  for (i in 1:length(my.files)){
    my.file <- read.table(my.files[i], header=TRUE)
    ifelse(i == 1, full.file <- my.file, full.file <- cbind(full.file, my.file[,c("beta","se","p_score")]))
           ifelse(i == 1, names(full.file)[9] <- paste(my.names[i], "_beta", sep=""), names(full.file)[(ncol(full.file)-2)] <- paste(my.names[i], "_beta", sep=""))
           ifelse(i == 1, names(full.file)[10] <- paste(my.names[i], "_se", sep=""), names(full.file)[(ncol(full.file)-1)] <- paste(my.names[i], "_se", sep=""))
           ifelse(i == 1, names(full.file)[13] <- paste(my.names[i], "_pscore", sep=""), names(full.file)[(ncol(full.file))] <- paste(my.names[i], "_pscore", sep=""))
  }
full.file <- full.file[,-c(2,11,12)]

#sort dataframe rows in order of Chrom, then Pos
str(full.file)
#full.file$ps <- as.numeric(full.file$ps)
#full.file$chr <- as.numeric(full.file$chr)
full.file <- full.file[with(full.file, order(chr, ps)), ]

#Make plotting variables
full.file$Index = NA
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(full.file$chr)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    full.file[full.file$chr==i, ]$Index=full.file[full.file$chr==i, ]$ps
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(full.file,full.file$chr==i-1)$ps, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    full.file[full.file$chr==i, ]$Index=full.file[full.file$chr==i, ]$ps+lastbase
  }
}

pheno.bin <- full.file
 #now, PosPhenos = 1 if p < avg 2500th SNP 1000x thr (99%), 0 if p > thr
 # also calculate/ report in text for 99.9% Thr = 250

##check threshold here
gethr <- mythrs[mythrs$SNPnum==2500,]
#get an ordered threshold list for the phenotypes
mythrlist <- rbind(gethr[gethr$pheno=="LA0410",],gethr[gethr$pheno=="LA3475",],gethr[gethr$pheno=="LA4345",],gethr[gethr$pheno=="LA4355",],gethr[gethr$pheno=="LA0480",],gethr[gethr$pheno=="LA1547",],gethr[gethr$pheno=="LA1589",],gethr[gethr$pheno=="LA1684",],gethr[gethr$pheno=="LA2093",],gethr[gethr$pheno=="LA2176",],gethr[gethr$pheno=="LA2706",],gethr[gethr$pheno=="LA3008",])

#remove Domesticated for this part
pheno.bin <- pheno.bin[,c(1:7,8:19,29:53)]
  ##unclear which phenotypes are significant without permutation analysis -- need to wait on this step until determining permutation threshold.
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")
  for (i in c(1:12)){
    mybeta = 5 + (3*i)
    mypscore = 7 + (3*i)
  pheno.bin[,paste(colnames(pheno.bin[mybeta]),"bin", sep="")] <- ifelse(pheno.bin[,mypscore] < mythrlist[i,3], 1, 0)
  }

names(pheno.bin)
pheno.bin$SUMM <- rowSums(pheno.bin[,c(45:56)], na.rm=T)

## check file name here
#write.csv(pheno.bin, "D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut99Thr_kmat1.csv")
#pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut99Thr_kmat1.csv")
table(pheno.bin$SUMM)
  
#high overlap SNP list for annotation
HOSNP <- pheno.bin[pheno.bin$SUMM > 5,]
## check file name here
#write.csv(HOSNP, "D_08_results/12Plants_HiOverlapSNPs_trueMAF20_10NA_GEMMA_1kpermut99Thr_kmat1.csv")

#and top 1000 SNPs > Thr per genotype
#then wide format
top1ksnp <- pheno.bin
#*_betabin columns tells SNP p < Thr or not
table(top1ksnp$`9_LA3008_betabin`)
#SUMM column tells which SNPs have no sig phenos
#and number of sig phenos per SNP
#sig phenos only p < 0.01:
top1ksnp <- top1ksnp[top1ksnp$SUMM>0,]
#now only keep rows with top 1000 per phenotype

for (y in c(10,13,16,19,22,25,28,31,34,37,40,43)){
  #sort data frame with small to high p values for each phenotype
  top1ksnp <- top1ksnp[order(top1ksnp[,y]),]
  #then reassign beta to be snp rank 
  top1ksnp[,(y-2)] <- 1:nrow(top1ksnp)
}

names(top1ksnp)

#only beta <= 1000 are significant
top1ksnp$top1k <- ifelse(top1ksnp$`1_LA0410_beta` < 1001 | top1ksnp$`2_LA0480_beta` < 1001 | top1ksnp$`3_LA1547_beta` < 1001 | top1ksnp$`4_LA1589_beta` < 1001 | top1ksnp$`5_LA1684_beta` < 1001 | top1ksnp$`6_LA2093_beta` < 1001 | top1ksnp$`7_LA2176_beta` < 1001 | top1ksnp$`8_LA2706_beta` < 1001 | top1ksnp$`9_LA3008_beta` < 1001 | top1ksnp$`10_LA3475_beta` < 1001 | top1ksnp$`11_LA4345_beta` < 1001 | top1ksnp$`12_LA4355_beta` < 1001, "topset", "omit")
table(top1ksnp$top1k)

top1ksnp <- top1ksnp[top1ksnp$top1k == "topset",]
#write.csv(top1ksnp, "D_08_results/12Plants_top1kSNPs_MAF20_10NA_GEMMA_kmat1_99Thr.csv")

#------------------------------------------------------------------  
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut999Thr_kmat1.csv")

#get thresholds here 
mythrs <- read.csv("D_07_randOUTS/GEMMA_1krand_thresholds.csv")
#now, PosPhenos = 1 if p < avg 2500th SNP 1000x thr (99%), 0 if p > thr
# also calculate/ report in text for 99.9% Thr = 250
gethr <- mythrs[mythrs$SNPnum==250,] #too  many at 99% level --> using 99.9%

#find top SNPs to mark: pheno.bin$SUMM and LA2093 pscore
pheno.bin.top <- pheno.bin[pheno.bin$X6_LA2093_pscore < gethr[8,3], ]
library(plyr)
#pheno.bin.top <- head(arrange(pheno.bin.top, desc(X6_LA2093_beta)), n=100)
pheno.bin.top <- pheno.bin.top[pheno.bin.top$SUMM>6,] #too  many at SUMM > 6 --> using SUMM > 10
pheno.bin.top$Index

#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(myGEMMA$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#get top SNPs for lines
pheno.bin.top$Index

table(pheno.bin$SUMM)

#make plots 
  setwd("~/Projects/BcSolGWAS")
  jpeg("paper/plots/addGEMMA/FigS3b_sigphenos_ManhattanPlot.jpg", width=8, height=5, units='in', res=600)
  ggplot(pheno.bin, aes(x=Index, y=SUMM))+
    colScale+ #remove for rainbow plot
    theme_bw()+
    #    scale_x_continuous(breaks = ticks)+
    geom_point(aes(color = factor(chr)))+
    labs(list(y="Plant Genotypes per Significant SNP", x="Chromosome position"))+
    #nrow=8
    theme(legend.position="none")+
    theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
    guides(col = guide_legend(nrow = 8, title=element_blank()))+
    scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    scale_y_continuous(breaks= c(0,2,4,6,8,10,12), labels=c("0","2","4","6","8","10","12"))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    geom_vline(xintercept=1885377, lty=2)+
    geom_vline(xintercept=2061494, lty=2)+
    geom_vline(xintercept=16687404, lty=2)+
    geom_vline(xintercept=17187675, lty=2)+
    geom_vline(xintercept=20864412, lty=2)+
    geom_vline(xintercept=22288525, lty=2)+
    geom_vline(xintercept=22783105, lty=2)+
    geom_vline(xintercept=25359723, lty=2)+
    geom_vline(xintercept=31953496, lty=2)+
    geom_vline(xintercept=34388621, lty=2)+
    geom_vline(xintercept=36077967, lty=2)
  dev.off()
  
# #overlaps for 99% threshold:
#   geom_vline(xintercept=1447936, lty=2)+
#     geom_vline(xintercept=1885377, lty=2)+
#     geom_vline(xintercept=2209430, lty=2)+
#     geom_vline(xintercept=4790255, lty=2)+
#     geom_vline(xintercept=7302540, lty=2)+
#     geom_vline(xintercept=12546559, lty=2)+
#     geom_vline(xintercept=13007560, lty=2)+
#     geom_vline(xintercept=16687404, lty=2)+
#     geom_vline(xintercept=20864412, lty=2)+
#     geom_vline(xintercept=21549298, lty=2)+
#     geom_vline(xintercept=21848271, lty=2)+
#     geom_vline(xintercept=22288525, lty=2)+
#     geom_vline(xintercept=22783105, lty=2)+
#     geom_vline(xintercept=25601190, lty=2)+
#     geom_vline(xintercept=26577239, lty=2)+
#     geom_vline(xintercept=31063497, lty=2)+
#     geom_vline(xintercept=31953496, lty=2)+
#     geom_vline(xintercept=34388621, lty=2)+
#     geom_vline(xintercept=36077967, lty=2)

#add figure 5a / b for GEMMA (S3): SNP overlap across phenotypes
table(pheno.bin$SUMM)

##check output file name
jpeg("paper/plots/addGEMMA/S3a_topSNPOverlap_12Plants_GEMMA_999thr.jpg", width=8, height=5, units='in', res=600)
ggplot(pheno.bin, aes(pheno.bin$SUMM)) + 
  geom_bar()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_y_continuous(name= "Number of SNPs", limits = c(0,20000))+
  scale_y_continuous(name= "Number of SNPs", limits = c(0,4500))+
  scale_x_continuous(name= "Plant Genotypes per Candidate SNP", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12), limits = c(0, 12))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
dev.off()


#small inset plot
##check output filename
jpeg("paper/plots/addGEMMA/S3a_topSNPOverlap_12Plants_GEMMA_inset_999thr.jpg", width=4, height=3, units='in', res=600)
ggplot(pheno.bin, aes(pheno.bin$SUMM)) + 
  geom_bar()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_y_continuous(name= "", limits = c(0,400))+
  scale_y_continuous(name= "", limits = c(0,25))+
  scale_x_continuous(breaks=c(6,7,8,9,10,11,12),labels=c(6,7,8,9,10,11,12), limits = c(5, 13), name="")+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
dev.off()

setwd("~/Projects/BcSolGWAS/paper/plots/addGEMMA")
#can add this later if needed. calculation needs fixing.
myprobs <- read.csv("GEMMA_probabilities.csv")

for (i in 44:55){
  myg <- sum(pheno.bin[,i])
  print (i)
  print(myg)}