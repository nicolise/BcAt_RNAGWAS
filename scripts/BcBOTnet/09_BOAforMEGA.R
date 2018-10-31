#Nicole E Soltis
#--------------------------------------------------

#need FASTA file input for MEGA
#want: beginning of C1 (5133) to 97414
#BcBoa17 goes up to 69449. 3 genes downstream is gene:Bcin01g00190. end position is 97414. 

#or, if looking within deletion ONLY:
#last base of deletion = 82614
#first base of deletion = first base of C1 = 4029

#approach 2: read in BOA data, convert to fasta format in R
#have been using binary PED -- lost SNP state info.
#need to use: 

install.packages("seqRFLP")
library("seqRFLP")
#Your sequences need to be in a data frame with sequence headers in column 1 and sequences in column 2 [doesn't matter if it's nucleotide or amino acid]
names <- c("seq1","seq2","seq3","seq4")
sequences<-c("EPTFYQNPQFSVTLDKR","SLLEDPCYIGLR","YEVLESVQNYDTGVAK","VLGALDLGDNYR")
df <- data.frame(names,sequences)
#Then convert the data frame to .fasta format using the function: 'dataframe2fas'
df.fasta = dataframe2fas(df, file="df.fasta")


#---------------------------------------------------------------------
#failed approach 1: read in Chr 1 fasta, try to extract nucleotides of interest
#BUT this does not include SNP state info for individual isolates. Dropping this approach.

setwd("~/Projects/BcGenome/data/ensembl/B05.10/fasta/singleChr")

library("Biostrings")
library("seqinr")

coi.fa <- read.fasta(file = textConnection("Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa/gpfs/nobackup/ensembl/amonida/rel91/eg38/fungi/fasta/botrytis_cinerea/dna/Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa"), as.string = T)

myFA <- read.fasta(file = "Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa/gpfs/nobackup/ensembl/amonida/rel91/eg38/fungi/fasta/botrytis_cinerea/dna/Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa", as.string=FALSE, seqonly=FALSE)
myFAdf <- as.data.frame(myFA)

myFA2 = readDNAStringSet("Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa/gpfs/nobackup/ensembl/amonida/rel91/eg38/fungi/fasta/botrytis_cinerea/dna/Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa")
seq_name = names(myFA2)
sequence = paste(myFA2)
myFAdf <- data.frame(seq_name, sequence)


library("Biostrings")
fasta2dataframe=function(fastaFile){
  s = readDNAStringSet(fastaFile)
  RefSeqID = names(s)
  RefSeqID = sub(" .*", "", RefSeqID) 
  #erase all characters after the first space: regular expression matches a space followed by any sequence of characters and sub replaces that with a string having zero  characters 
  
  for (i in 1:length(s)){
    seq[i]=toString(s[i])
  }
  
  RefSeqID_seq=data.frame(RefSeqID,seq)
  return(RefSeqID_seq)
}
mydf2 = fasta2dataframe("Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa/gpfs/nobackup/ensembl/amonida/rel91/eg38/fungi/fasta/botrytis_cinerea/dna/Botrytis_cinerea.ASM83294v1.dna.chromosome.1.fa")
