rm -r ~/hg-data/ ~/hg-tmp/
higlass-manage start -n Davos -p 8888

for c in $(ls ../../Data/Processed/2019-06-18_PCa-LowC-sequencing/Contacts/); do
    echo $c;
    higlass-manage ingest --hg-name Davos --filetype cooler --datatype matrix --assembly hg38 --project-name Davos $c
done

for b in $(ls ../2020-01-15_TAD-aggregation/resolved-TADs/PCa*.bedGraph); do
    higlass-manage ingest --hg-name Davos --filetype bedfile --datatype bedlike --assembly hg38 $b
done