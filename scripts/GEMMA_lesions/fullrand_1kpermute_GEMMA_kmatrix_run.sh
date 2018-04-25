#!/bin/sh

#script is run in B05_GEMMA_les/
#files are structured as B05_GEMMA_les/D_05_bigrand/rand1k_1 ... through 1000)
#only run first 4 lesion phenotypes for permutations: Col0, coi1, npr1, pad3

for i in {1..1000}
do
  mybpath='D_05_bigrand/rand1k_'"$i"'/binMAF20NA10_rand' 
  myoutpath='rand1k_'"$i"
  mkdir -p ./output/$myoutpath;

    for j in {1..4}
    do
        echo "Looping ... phenotype $j"
        ~/gemma/bin/gemma -bfile $mybpath -k D_03_kmat/binMAF20NA10_PLINK_randtest_kmatrix1.cXX.txt -n $j -miss 0.1 -maf 0.2 -lm 4 -o $myoutpath"/pheno${j}"
    done
done
