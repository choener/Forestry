#!/bin/bash

RNAfold="/scratch/choener/Software/ViennaRNA-2.3.5/bin/RNAfold"

# which sequence lengths to use
for i in `seq 25 25 200`
do
  printf "%4d " $i
  zcat examples/0300.input.gz | cut -c-$i | \
    $RNAfold | \
    time -f "%U" "$@" \
    > /dev/null \
    2>tmp.tmp
    SEC=`cat tmp.tmp`
    printf "%0.3f\n" $(bc -l <<< "$SEC / 50")
done
