#Nicole E Soltis
#05/29/18
#testing R meta-analysis tools from GWAS studies
#------------------------------------------------------------------------------
#sample data from BcSolGWAS - can use 12phenos

library(catmap) #may only work for case/ control or family-based data
library(rmeta) #appears to do only treatment/ control
library(metafor) #could work... has it been used for GWA?

#first must calculate effect sizes
?escalc()
#assumes meta analysis types: across-group comparisons (~anova?), 2-variable correlation (~linear regression?) or other. Unsure how to code GWA...

#want: forest plot for GWA
?forest()

#summarizes "true effect" based on observed effect size across each of i studies
#may need to focus this on just putative trans hotspots?? top hits per gene outside of cis region?? idk
