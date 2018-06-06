#Nicole E Soltis adapted from Andrey A. Shabalin
#code adapted from http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
#05/22/18

#---------------------------------------------------------------------------------
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

#file format: genotype, expression, covariates, gene location, SNP location
#first 3: long format. each row is one observation. column orders must match. check as.numeric()
#genotype: cols are (SNP) id, then individual 1...n. SNP states are 0, 1, 2 (1 is het)
#expression: (gene) id, then individual 1...n with expression levels
#covariates: id, then individual 1...n (e.g. sex, age...)
#gene loc: cols are geneid (match expression ids), chr (format "chr1"), s1, s2 (start and stop position ??)
#SNP loc: snp (match genotype ids), chr (format "chr1"), pos

#Matrix eQTL is designed to handle large genotype and expression data sets. They are loaded using SlicedData classes which store the data in slices of 1000 rows (default size). The analysis is then performed for each pair of slices of genotype and expression data sets.
#advised to use fast BLAS: best done in Linux!

#Linux version of R: consult the R Installation and Administration manualexternal link and a relevant blog postexternal link by Allan Engelhardt. 
#https://cran.r-project.org/doc/manuals/R-admin.html#BLAS
#https://www.r-bloggers.com/faster-r-through-better-blas/ 

#error types, FDRs, etc: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/features.html 
