#021518
#Nicole E Soltis
#Compare_phenotype_methods.R

#test different scaling approaches for low lsmeans transcripts

#-----------------------------------------------------------------
rm(list=ls())

#load data: individual plant genotypes. Directly from Vivian.
setwd("~/Documents/GitRepos/BcAt_RNAGWAS/data")
setwd("~/Projects/BcAt_RNAGWAS/data")

#load data: lsmeans
#lsm.for.GWAS comes from allreads/01_PrepPhenos_lsmeans.R
lsm.for.GWAS <- read.csv("allreadsGWAS/T4/01_prepFiles/lsmeans_allreads.csv")

#keep only 10ish transcripts with 10% or more isolates with lsmeans < -20
#10% of 96 is 10 reads < -20
iso.names <- lsm.for.GWAS[,2]
all.my.reads <- lsm.for.GWAS[,3:9269]
my.min <- all.my.reads[1,]
for (i in (1:9267)){
  for (j in (1:10)){
    my.min[j,i] <- sort(all.my.reads[,i],partial=j)[j]   
  }
}
my.min.lowest <- my.min[,my.min[5,]<(-26)] #I'll just use this low group of 12 transcripts

#minimum:Bcin11g00300.1
my.names <- "Bcin11g00300.1"
my.low.reads <- as.data.frame(lsm.for.GWAS[,names(lsm.for.GWAS) %in% my.names])
my.low.reads$Isolate <- lsm.for.GWAS[,2]
names(my.low.reads)[1] <- "Bcin11g00300.1"
my.low.reads$Bcin11g00300.1 <- as.numeric(my.low.reads$Bcin11g00300.1)

# my.low.reads.cutoff <- my.low.reads
# for (y in (2:13)){
#   my.low.reads.cutoff[,y] <- ifelse(my.low.reads.cutoff[,y] < (-5), (-5), my.low.reads.cutoff[,y])
#   names(my.low.reads.cutoff)[y] <- paste(names(my.low.reads.cutoff)[y], "cutoff", sep=".")
# }

my.low.reads.z <- my.low.reads
for (z in c(1)){
  my.low.reads.z[,z] <- scale(my.low.reads.z[,z], center = TRUE, scale = TRUE)
  names(my.low.reads.z)[z] <- paste(names(my.low.reads.z)[z], "zscale", sep=".")
}

IsoNames <- read.csv("IsolateKey_Vivian.csv")
MyReads <- read.csv("Vivian_Bc/result.lsm.csv")
MyReads <- MyReads[,-c(1)]

#attach Isolate Names
names(IsoNames)[1] <- "Isolate"
IsoNames <- IsoNames[,c(1,3)]
MyReads <- merge(MyReads, IsoNames, by="Isolate")
MyReads <- MyReads[,c(9270,2,3:9269)]
names(MyReads)[1] <- "Isolate"
#only keep target genotype Bcin11g00300.1
MyReads <- MyReads[,c("Isolate","HostGenotype","Bcin11g00300.1")]
names(MyReads)[2]<- "Phenotype"

#each plant x transcript is a unique phenotype
#split by host
splitlist <- split(MyReads, MyReads$Phenotype)
MyReads.coi1 <- splitlist[[1]]
names(MyReads.coi1)[3] <- "Mean Count"
MyReads.col0 <- splitlist[[2]]
names(MyReads.col0)[3] <- "Mean Count"
MyReads.npr1 <- splitlist[[3]]
names(MyReads.npr1)[3] <- "Mean Count"


#merge
my.low.reads$Phenotype <- "lsmeans"
names(my.low.reads)[1] <- "Mean Count"
my.low.reads <- my.low.reads[,c(2,3,1)]
my.low.reads.z$Phenotype <- "lsmeans_z"
names(my.low.reads.z)[1] <- "Mean Count"
my.low.reads.z <- my.low.reads.z[,c(2,3,1)]

my.plotdat <- rbind(my.low.reads, my.low.reads.z, MyReads.col0, MyReads.coi1, MyReads.npr1)
names(my.plotdat)[3] <- "Mean_Count"

setwd("~/Projects/BcAt_RNAGWAS/plots")
jpeg("testMethods/comparelsmeans_violin_1gene.jpg", width=8, height=5, units='in', res=600)
x$name <- factor(x$name, levels = x$name[order(x$val)])
my.plotdat$Phenotype <- factor(my.plotdat$Phenotype, levels = c("col.0","coi.1","npr.1","lsmeans","lsmeans_z"))
p <- ggplot(my.plotdat, aes(factor(Phenotype), Mean_Count))
p + geom_violin(fill="slateblue3")+
  theme_bw()+
  labs(list(y="Mean Read Count", x=""))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
dev.off()
