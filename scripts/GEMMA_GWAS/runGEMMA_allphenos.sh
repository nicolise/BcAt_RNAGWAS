#!/bin/sh
for i in {1..10}
#in reality, 6:9272
do
  #j = $(($i - 5))
  #echo "Looping ... number $j"
  echo "Looping ... number $i"
  ~/gemma/bin/gemma -bfile 02_GEMMA/binMAF20NA10 -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_PLINKv3"_${i}"
done
