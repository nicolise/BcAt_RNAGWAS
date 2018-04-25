#!/bin/sh
  echo "only running once for default phenotype: previous runs on DWS told us the k-matrices for all phenos are same!"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile D_02_randGEMMA/binMAF20NA10 -gk 1 -o binMAF20NA10_PLINK_randtest_kmatrix1
  echo "k matrix 2"
  ~/gemma/bin/gemma -bfile D_02_randGEMMA/binMAF20NA10 -gk 2 -o binMAF20NA10_PLINK_randtest_kmatrix2
done
