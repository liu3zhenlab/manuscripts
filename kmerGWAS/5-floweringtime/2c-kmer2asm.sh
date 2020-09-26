#!/bin/bash
cap3=/homes/liu3zhen/software/cap3/CAP3/cap3
B73Ref4=/homes/liu3zhen/references/B73Ref4/blast+/B73Ref4
A188Ref1=/homes/liu3zhen/references/A188Ref1/genome/blast+/A188Ref1
for km in 1o-clust*.fasta; do
	$cap3 $km -i 21 -j 31 -s 251 -o 16 -p 90 1>2o-cap3.log 2>&1
done

