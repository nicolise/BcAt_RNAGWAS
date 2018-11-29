#Nicole E Soltis
#11/05/18

#-----------------------------------------------------------------------------

rm(list=ls())
#plot/ table: How often is max real p < max permut p?
setwd("~/Projects/BcAt_RNAGWAS/data/GEMMA_eachAt_Bc/")
randdat <- read.csv("06_GEMMAsumm_RAND/TranscNames/col0_GEMMA_top1SNPsample_genes.csv")
#need annotated gene names ("phenotype") to match random to real data

mydat01 <- read.csv("06_GEMMAsumm/GeneNames/col0_GEMMA_top1SNPsample.csv")

randmerge <- randdat[,c("Gene","chr","ps", "p_score","randrun")]
names(randmerge) <- c("Gene","rand_chr","rand_pos", "rand_p","rand_run")
hotspt <- merge(mydat01, randmerge, by="Gene")
hotspt$DvR <- hotspt$rand_p - hotspt$p_score
hist(hotspt$DvR)
hotspt$mygroup <- ifelse(hotspt$DvR > 0, "NonSig", "Sig")
#add indexing now
mydat_plot <- hotspt
mydat_plot <- mydat_plot[order(mydat_plot$chr, mydat_plot$ps),]
mydat_plot$Index.s = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(mydat_plot$chr)) {
  print(i)
  if (i==1) {
    mydat_plot[mydat_plot$chr==i, ]$Index.s=mydat_plot[mydat_plot$chr==i, ]$ps
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==i-1)$ps, 1)
    mydat_plot[mydat_plot$chr==i, ]$Index.s=mydat_plot[mydat_plot$chr==i, ]$ps+lastbase
  }
}

#--------------------------------------------------------------------
#mini manhattan: focus in on chr 9
#plot max rand on Chr 9 only. Label Mb X axis. 
setwd("~/Projects/BcAt_RNAGWAS")
mydat_c9 <- mydat_plot[mydat_plot$chr==9,]

#--------------------------------------------------------------------------
#could also plot this as distance -- set all positive values to zero, then plot -1*x
mydat_plot$DvR <- ifelse(mydat_plot$DvR > 0, 0, mydat_plot$DvR)

#-----------------------------------------------------------------------
#plot with DF split between groups (sig, nonsig)
#split df by groups - get mydat_plot from chunk above
split.hot <- split(mydat_plot, mydat_plot$mygroup)
assign("hotspt.ns", split.hot[[1]])
assign("hotspt.sig", split.hot[[2]])

## repeat this once each for .sig, .ns
#mydat_plot <- hotspt.sig
mydat_plot <- hotspt.sig

#summarize within each SNP - # of transcript hits
mydat_summ <- mydat_plot[,c("chr","ps","p_score","Gene","Index.s")]
mydat_summ_ngene <- aggregate(Gene ~ Index.s, data = mydat_summ, FUN = function(x){NROW(x)})
#now add SNP data back on, matching by Index.s
mydat_summ_ngene <- merge(mydat_summ_ngene, mydat_summ[,c("chr","ps","Index.s")], by="Index.s")
#remove duplicate rows
mydat_summ_ngene <- unique(mydat_summ_ngene)
mydat_plot <- mydat_summ_ngene

#-------------------------------------------------
#troubleshooting here
#check chromosome 12 high-p SNPs-- no overlap?
mydat01_c12 <- mydat01[mydat01$chr==12,]
mydat_plot_c12 <- mydat_plot[mydat_plot$chr==12,]
hotspt.sig_c12 <- hotspt.sig[hotspt.sig$chr==12,]
hotspt.ns_c12 <- hotspt.ns[hotspt.ns$chr==12,]
