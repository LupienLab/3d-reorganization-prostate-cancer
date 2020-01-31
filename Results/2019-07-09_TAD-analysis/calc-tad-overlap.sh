#!/bin/bash

SAMPLE_IDS=($(eval echo "PCa{13266,13848,14121,19121,3023,33173,40507,51852,53687,56413,57054,57294,58215}"))
mkdir -p TAD-comparisons/intersected/

for ((w=3; w <= 30; w++));
do
    echo "${w}"
    for ((i=0; i < 13; i++));
    do
        s1=${SAMPLE_IDS[$i]}
        echo "\t${s1}"
        for ((j=i + 1; j < 13; j++));
        do
            s2=${SAMPLE_IDS[$j]}
            echo "\t\t${s2}"
            # only take TADs whos window size parameter corresponds to the current w being looped over
            bedtools intersect -f 0.6 -r -a ../2020-01-15_TAD-aggregation/resolved-TADs/separated-TADs/${s1}.40000bp.w_${w}.domains.bed -b ../2020-01-15_TAD-aggregation/resolved-TADs/separated-TADs/${s2}.40000bp.w_${w}.domains.bed > TAD-comparisons/intersected/w_${w}.${s1}.${s2}.intersected.bed
        done
    done
done
