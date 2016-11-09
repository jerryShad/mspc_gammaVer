# MSPC - Bioconductor Package for Multiple Sample Peak Calling

The analysis of ChIP-seq samples outputs a number of enriched regions, each indicating a protein-DNA interaction or a specific chromatin modification. Enriched regions (commonly known as "peaks") are called when the read distribution is significantly different from the background and its corresponding significance measure (p-value) is below a user-defined threshold.

When replicate samples are analysed, overlapping enriched regions are expected. This repeated evidence can therefore be used to locally lower the minimum significance required to accept a peak. Here, we propose a method for joint analysis of weak peaks.

Given a set of peaks from (biological or technical) replicates, the method combines the p-values of overlapping enriched regions: users can choose a threshold on the combined significance of overlapping peaks and set a minimum number of replicates where the overlapping peaks should be present. The method allows the "rescue" of weak peaks occuring in more than one replicate and outputs a new set of enriched regions for each replicate.
