
module load vcftools


# download fastPhase

wget http://scheet.org/code/Linuxfp.tar.gz
tar -zxvf Linuxfp.tar.gz 


## get conversion script

wget https://github.com/lstevison/vcf-conversion-tools/blob/master/vcf2fastPHASE.pl

chmod 755 vcf2fastPHASE.pl



## ONE STEP NOT IN HERE- Filter vcf to be region of interest! Easy to do with vcftools


## filter biallelic

vcftools --vcf test_te.recode.vcf --min-alleles 2 --max-alleles 2 --recode --out test_te.biallelic.vcf


## convert to fastphase format

## arg 1 is input vcf

## arg2 and 3 are output files, geno one is what you give to fastphase

## arg4 is number individuals (this is output when run vcftools above)

perl vcf2fastPHASE.pl test_te.biallelic.vcf.recode.vcf test_te.fastPHASE.genos test_te.fastPHASE.positions 1135


## run fastphase

# the Pzp option outputs the expected total cluster membership for each snp

./fastPHASE -T1 -Pzp test_te.fastPHASE.genos 
