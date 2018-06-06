#!/bin/sh
for i in {1..2}
do
  echo "Looping ... phenotype $i"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile D_02_randGEMMA/binMAF20NA10 -k D_03_kmat/binMAF20NA10_PLINK_randtest_kmatrix1.cXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_fullrand_kmat1"_pheno${i}"
  echo "k matrix 2"
  ~/gemma/bin/gemma -bfile D_02_randGEMMA/binMAF20NA10 -k D_03_kmat/binMAF20NA10_PLINK_randtest_kmatrix2.sXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_fullrand_kmat2"_pheno${i}"
  echo "no k matrix"
  ~/gemma/bin/gemma -bfile D_02_randGEMMA/binMAF20NA10 -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_fullrand_kmat0"_pheno${i}"
done
