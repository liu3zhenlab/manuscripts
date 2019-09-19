### Effector Gene Reshuffling Involves Dispensable Mini-chromosomes in the Wheat Blast Fungus
In-house scripts used in the analyses are posted here.

1. dotplot R script to compare two sets of sequences
first, nucmer was used to generate alignments
```
#!/bin/bash

query=B71Ref1.fas
mg8=Magnaporthe_oryzae.MG8.dna.toplevel.fa
nucmer --maxmatch --nosimplify $mg8 $query -l 1000 -p 3o-B71Ref1toMG8
```
second, the function dotplot.two.sets.seqs.R is loaded and run
```
# R scripts
source("dotplot.two.sets.seqs.R")
deltafile <- "3o-B71Ref1toMG8.delta"
nucmer.twosets.plot(datapath=".", datafile=deltafile, lend.turnoff=F,
  aln.min.size=10000, aln.min.identity=95,
  tableout=T, xlab = "MG8", ylab = "B71Ref1",
  main = "B71Ref1 vs.MG8",
  line.width.factor=2,
  xlabel.rm = "supercont", ylabel.rm = "chr",
  pdfout=T, outpath=".", pdf.width = 4.9, pdf.height = 5,
  tableoutfile="B71Ref1toMG8.nucmer.txt",
  imageoutfile="FigDotPlot.B71Ref1toMG8.nucmer.pdf")
```
 
2. R script [DESeq.single.trt.R](DESeq.single.trt.R) is to used to detect differential expression between two groups.
