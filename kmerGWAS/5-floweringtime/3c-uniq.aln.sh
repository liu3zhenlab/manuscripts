#!/bin/bash
A188Ref1=~/references/A188Ref1/genome/bowtie/A188Ref1
B73Ref4=~/references/B73Ref4/bowtie/B73Ref4
for kfasta in 1o*fasta.cap.contigs; do
	out=`echo $kfasta | sed 's/^1o-/3o-/g' | sed 's/.fasta.cap.contigs//g'`
	bowtie -k 1 -m 1 --best --strata -v 2 -f --sam-nohead -B 1 $A188Ref1 $kfasta > ${out}.A188uniq.sam
	bowtie -k 1 -m 1 --best --strata -v 2 -f --sam-nohead -B 1 $B73Ref4 $kfasta > ${out}.B73uniq.sam
done
