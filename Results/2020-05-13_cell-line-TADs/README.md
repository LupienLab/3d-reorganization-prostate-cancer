# Summary

This folder contains the TADs called by TopDom for various cell lines.

For many of these cell lines, there are few read pairs associated with `chrY`, which elads to errors in calling TADs on that chromosome.
To avoid modifying the `topdom.R` script, I have opted to ignore `chrY` TAD calls for the cell lines.

I will have to deal with this when comparing to our primary samples, since the code includes `chrY`, but that is a problem for the future.
