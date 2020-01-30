#!/bin/bash

SAMPLE_IDS=($(eval echo "PCa{13266,13848,14121,19121,3023,33173,40507,51852,53687,56413,57054,57294,58215}"))

for ((w=3; w < 30; w++));
do
    echo "${w}"
    for ((i=0; i < 13; i++));
    do
        s1=${SAMPLE_IDS[$i]}
        echo -e "\t${s1}"
        for ((j=i + 1; j < 13; j++));
        do
            s2=${SAMPLE_IDS[$j]}
            echo -e "\t\t${s2}"
            # only take TADs whos window size parameter corresponds to the current w being looped over
            bedtools intersect \
                -f 0.6 \
                -r \
                -a <(awk -v w=${w} '{if ($6 == w) print}' ../2020-01-15_TAD-aggregation/resolved-TADs/${s1}.40000bp.aggregated-domains.sorted.bedGraph) \
                -b <(awk -v w=${w} '{if ($6 == w) print}' ../2020-01-15_TAD-aggregation/resolved-TADs/${s2}.40000bp.aggregated-domains.sorted.bedGraph) \
                > TAD-comparisons/w_${w}.${s1}.${s2}.intersected.bed
            cp TAD-comparisons/w_${w}.${s1}.${s2}.intersected.bed TAD-comparisons/w_${w}.${s2}.${s1}.intersected.bed
        done
    done
done
