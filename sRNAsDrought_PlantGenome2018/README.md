## sRNAsDrought_PlantGenome2018
keys scripts used in the manuscript of maize small RNA transcriptomic profiles under drought. The manuscript was published in the Plant Genome (doi: 10.3835/plantgenome2018.08.0058).

This directory includes key scripts for analyzing a time-series small RNA (sRNA) transcriptomic data to identify and cluster drought-responsive sRNAs.

1. DESeq2 analysis to identify drought-responsive sRNAs  
Rscript sRNA.DESeq2.interaction.R  
**Note.** sRNA.counts.txt was used as the input data in the analysis.


2. mclust to cluster drought-responsive sRNAs  
Rscript sRNA.mclust.R  
**Note.** Normalized counts of drought-responsive sRNAs in the file of sig.sRNA.normalized.counts.txt were used as the input data.
