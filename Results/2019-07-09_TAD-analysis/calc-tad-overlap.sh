for w in {3..30}; do
    echo "${w}"
    for s1 in PCa{3023,13266,13848,14121,19121,33173,40507,51852,53687,56413,57054,57294,58215}; do
        echo -e "\t${s1}"
        for s2 in PCa{3023,13266,13848,14121,19121,33173,40507,51852,53687,56413,57054,57294,58215}; do
            echo -e "\t\t${s2}"
            bedtools intersect -f 0.6 -r -a TADs/w_${w}/${s1}.40000bp.domains.bed -b TADs/w_${w}/${s2}.40000bp.domains.bed > TADs/Comparisons/w_${w}.${s1}.${s2}.intersected.bed
        done
    done
done
