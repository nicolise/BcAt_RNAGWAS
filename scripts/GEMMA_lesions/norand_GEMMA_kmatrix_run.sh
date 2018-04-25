#!/bin/sh
for i in {1..6}
do
  echo "Looping ... phenotype $i"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile D_02_randGEMMA/binMAF20NA10_rand -k D_03_kmat/binMAF20NA10_PLINK_randtest_kmatrix1.cXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_fullrand_kmat1"_pheno${i}"
done
