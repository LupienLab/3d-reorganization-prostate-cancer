# Summary

In `2019-07-24_breakfinder/`, with default parameters we found that the positive control for our data (the TMPRSS2-ERG fusion) was only detected in one sample, PCa51852.
It had the second highest number of reads of the T2E samples, only behind PCa57054, which has an increased expression of ERG, but not the T2E fusion.

From the contact matrices, it appears that it should be detected in all samples, but wasn't.
To rule out sequencing depth as a contributing factor in the detection of the fusion, we're subsampling the PCa51852 BAM to the same depth as the shallowest sample (453.8 M unique read pairs in PCa13848), and running breakfinder again at a variety of thresholds to detect the fusion.

## Results

After downsampling the BAM for PCa51852 5 separate times (5 different seeds, all taking 58%), each run of Breakfinder with default parameters has called the T2E fusion.

## Conclusions

It appears that sequencing depth is not an issue for detection, for some reason, despite not being detected in the other samples.
It is unclear why this fusion is not being called in the remaining samples, then.
