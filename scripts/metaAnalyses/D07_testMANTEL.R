#Nicole E Soltis
#R script to prep data for MANTEL post-hoc p-value meta-analysis
#06/06/18
#----------------------------------------------------------------
#see MANTEL.pl for file instructions
#maybe not ideal: 8n + 1 columns in dataframe for n phenotypes, x rows for x SNPs
###################################################################
#
# Required input:
#
# 1) A plain-text file with association analysis results for all studies combined into a single file
#    (one line per SNP).
#
#    Expected file format (where the number indicates the column):
#
#    SNP CHR POS BETA SE PVALUE CA A2 CAF RATIO
#      1   2   3    4  5      6  7  8   9    10
#
#    where CHR is the chromosome number, POS is the chromosomal position of the SNP,
#    BETA is the computed estimate parameter of the effect of the given SNP,
#    SE is the standard error around that BETA estimate,
#    CA represents the coded allele that the BETA and SE are referring to, A2 is the other allele,
#    CAF refers to the allele frequency of the coded allele (CA), and
#    RATIO is the ratio of the observed variance of the dosage to the expected (binomial) variance.
#
#    Columns 2-10 are repeated (on the same line) for every additional GWAS that is part of the
#    meta-analysis.
#
#    The RATIO is used to correct the weight of the contribution of each individual study depending
#    on the imputation quality of the SNP.  This is only necessary when imputation was used (set to
#    1 if all SNPs are genotyped experimentally).   See de Bakker et al., Human Molecular Genetics,
#    2008 for more background information on this topic.
#
# 2) A plain-text file that contains study-specific parameters
#
#    Example file format:
#
#    FHS    1.023     7650      1
#    CHS    1.040      854     17.2755
#    ERGO   1.034     4606     18.1075
#
#    Column 1 : contains an alphanumeric name to identify the studies listed in the file specified
#               above
#    Column 2 : lists the genomic inflation factor (lambda) -- used for adjusting the SE on the fly
#               (set to 1 if the association results are already adjusted)
#    Column 3 : lists the sample size of each study -- used for meta-analysis based on sample size-
#               weighted z-scores
#    Column 4 : lists the correction factor to standardize the BETA and SE estimates across all
#               studies to ensure the scale and units of the BETA and SE are identical (the BETA
#               and SE are divided by the given factor)
#
#    Note that the order of the studies in this file is important -- it reflects the order of the
#    association results as they appear in the data file specified above.
#
# 3) A PLINK BIM file -- for SNP positions
#
#    Expected file format (note that a PLINK BIM file has NO header line!):
#
#    CHR  SNP     MORGAN     POS  A    B
#      1    2          3       4  5    6
#
# 4) A PLINK generated file (--freq) for HapMap which is used as the reference to resolve
#    ambiguities in allele coding.
#
#    Expected file format:
#
#    CHR          SNP   A1   A2          MAF  NCHROBS
#      1            2    3    4            5        6
#
# 5) Gene annotations
#
#    Expected file format:
#
#    CHR START STOP GENE_SYMBOL
#
# 6) Maximal distance to a gene (in kilobase units) -- this is used in the final output as for
#    every SNP genes within the specified distance are listed.  A reasonable choice is 100-200
#    kilobases.

#------------------------------------------------------
#files needed for run_mantel.csh: see linux MANTEL/NES_example

#snps: 1 column, no header
  #c1 = unique snp names, must match c2 of freq

#dbsnp: 8 columns, tab separated, headers
  #c1 = chrom (chr1 format)
  #c2 = chromStart
  #c3 = chromEnd
  #c4 = SNP (unique name or + if no data)
  #c5 = strand (+ or -)
  #c6 = observed (A/G or -/C format etc) 
  #c7 = class (can call all = single)
  #c8 = func (can call all = unknown)

#freq: 6 columns, tab separated, headers
  #c1 = CHR (numeric)
  #c2 = SNP (unique name, must match dbsnp c4)
  #c3 = A1 (SNP state allele 1, ACTG0)
  #c4 = MAF (numeric)
  #c5 = NCHROBS (what does this even mean? numeric, 114:222 ish)

#genes: 4 columns, space separated, no headers
  #c1 = chromosome, must match freq c1
  #c2 = snp start, c2 = snp end, c4 = gene name