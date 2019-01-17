#10_drawHotspots.R
#Nicole E Soltis
#01/17/19

#---------------------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS")
#Index MidGene to space out chromosomes
mydat_plot <- read.csv("data/GEMMA_eachAt_Bc/06_GEMMAsumm/GeneNames/SNPannot/col0_GEMMA_top1SNPsample_genes_indexed.csv")
mydat_index <- mydat_plot[,c("chr.snp","ps")]
mydat_index <- mydat_index %>%
  group_by(chr.snp) %>%
  summarize(midgene=max(ps, na.rm=TRUE))
mydat_index <- as.data.frame(mydat_index)
#hand copied this over
#--------------------------------------------------------------------------------
MyHots <- read.csv("data/GEMMA_eachAt_Bc/07_TopSNPs/BcAt_permut/BcAt_topHots_Gene_funcannot_index.csv")

#add conditional hotspot size variable
MyHots$NumGenes <- ifelse(MyHots$HotspotCategory=="Bc", MyHots$Gene.B,MyHots$Gene.A)
MyHots$HotspotCategory <- droplevels(MyHots$HotspotCategory)

#now add indexing
#Make plotting variables for snp
mydat_plot <- MyHots[order(MyHots$chr, MyHots$midgene),]
mydat_plot$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (j in unique(mydat_plot$chr)) {
  print(j)
  if (j==1) {
    mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$midgene
  }	else {
    lastbase=lastbase+max(subset(mydat_plot,mydat_plot$chr==j-1)$midgene, 1)
    mydat_plot[mydat_plot$chr==j, ]$Index=mydat_plot[mydat_plot$chr==j, ]$midgene+lastbase
  }
}


library(ggplot2)
myColors <- c("darkgreen", "navyblue")
names(myColors) <- levels(MyHots$HotspotCategory)
colScale <- scale_colour_manual(name = "Trait",values = myColors)

setwd("~/Projects/BcAt_RNAGWAS")
jpeg("plots/paper/BcAt_HotSpots_LinearPlot.jpg", width=10, height=2, units='in', res=600)
print(
  ggplot(mydat_plot, aes(x=Index, y=(0), size=NumGenes))+
    theme_bw()+
    colScale+
    geom_point(aes(color=factor(HotspotCategory),alpha=0.1))+
    labs(list(y=NULL, title=NULL))+
    theme(legend.position="none")+
    scale_size_continuous(range = c(2, 20))+
    scale_x_continuous(name="Chromosome", breaks = c(2029725, 5715883, 9002014, 11775203, 14410595, 17176482, 19845645, 22470978, 25004941, 27457400, 29808907, 32126298, 34406278, 36587755, 38708818, 40640197, 41655662, 41838837), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16", "17","18"))
)
dev.off()