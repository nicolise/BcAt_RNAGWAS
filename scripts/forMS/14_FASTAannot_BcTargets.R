#Nicole E Soltis

#get FASTA sequences for Boytrytis hotspot targets-- look for secretion enrichment etc. 
#--------------------------------------------------------------------------------
#supplemental data set 1, target gene column has gene names
rm(list=ls())
setwd("~/Projects/BcAt_RNAGWAS/")
mydat <- read.csv("paper/plots/SDS1.csv")

#having trouble getting the right fasta for B05.10 genes. Can do individual looku pon fungi.ensembl.org/Botrytis_cinerea/Gene and then select sequence --> download sequence BUT THAT TAKES A MILLION YEARS

#alternative: paste a chunk of gene names into http://botbioger.versailles.inra.fr/botportal/cgi-bin/display.py?displaytype=genequery and at the upper left corner of the chart, select "quick access to all ORF sequences". Then run in http://www.cbs.dtu.dk/services/SignalP/

#also, how did I do GO analysis for previous Bc genetics?

#-----------------------------------------------------------------------------------------------------
#read in fasta file
library("Biostrings")
setwd("~/Projects/BcSolGWAS")
FullFASTA <- readDNAStringSet("data/BcGenome/WGS/botrytis_cinerea__t4__1_genes.fasta")
tryFASTA <- readDNAStringSet("data/BcGenome/WGS/87_iso_no_organic.fasta.tar/87_iso_no_organic.fasta")
setwd("~/Projects/BcGenome")
tryFASTA <- readDNAStringSet("data/ensembl/B05.10/fasta/Bc_fullgenome_chr.fa")
FullFASTA <- tryFASTA
#turn it into a dataframe to make it easy to work with
seq_name <- names(FullFASTA)
sequence <- paste(FullFASTA)
FASTAdf <- data.frame(seq_name, sequence)
#now only take characters before | pipe into a new variable to match with SNPs "Gene"
names(mydat)
names(mydat)[2] <- "HotspotGene"
names(mydat)[3] <- "TargetGene"

x <- "BcT4_1 | Botrytis cinerea (T4) BcT4_1 (1968 nt)"
x <- gsub("\\|.*", "", x)
x <- gsub(" ", "", x, fixed = TRUE)
names(FASTAdf)
FASTAdf$Gene <- FASTAdf$seq_name
FASTAdf$Gene <- gsub("\\|.*", "", FASTAdf$Gene)
FASTAdf$Gene <- gsub(" ", "", FASTAdf$Gene, fixed = TRUE)
names(FASTAdf)
FASTAgetseqs <- FASTAdf[,c("Gene", "sequence")]

-log10(0.01)

#now add sequences onto SNP lists
names(DomestGenes)
DomestAnnot <- DomestGenes

#individual domestication phenotypes
stop("feed in selected phenotype here")
DomestAnnot.dom <- subset(DomestAnnot, DomestAnnot[ ,2] != 0) 
DomestAnnot.wild <- subset(DomestAnnot, DomestAnnot[,3] != 0)
DomestAnnot.sens <- subset(DomestAnnot, DomestAnnot[,4] != 0)

DomestAnnot <- DomestAnnot.dom

#all domestication phenotypes:
DomestAnnot <- as.data.frame(DomestAnnot[,c("geneID")])
DomestAnnot <- rename(DomestAnnot, "DomestAnnot[, c(\"geneID\")]" = "Gene")
DomestAnnot <- merge(DomestAnnot, FASTAgetseqs, by="Gene")

#now convert back to FASTA
names(DomestAnnot)
#install.packages("seqinr")
library(seqinr)
write.fasta(as.list(DomestAnnot$sequence),DomestAnnot$Gene,"data/GWAS_files/05_annotation/DomestGenes_domest_NA10.fasta")