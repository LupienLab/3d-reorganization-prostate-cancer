for w in {30,20,10,3}; do
    echo "${w}"
    for s1 in PCa{3023,13266,13848,14121,19121,33173,40507,51852,53687,56413,57054,57294,58215}; do
        echo -e "\t${s1}"
        for s2 in PCa{3023,13266,13848,14121,19121,33173,40507,51852,53687,56413,57054,57294,58215}; do
            echo -e "\t\t${s2}"
            for c in chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}; do
                echo -e "\t\t\t${c}"
                bedtools intersect -f 0.6 -r -a w_${w}/${s1}.40000bp.${c}.domains.bed -b w_${w}/${s2}.40000bp.${c}.domains.bed >> Comparisons/w_${w}.${s1}.${s2}.intersected.bed
            done
        done
    done
done
