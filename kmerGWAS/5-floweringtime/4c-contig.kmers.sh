#!/bin/bash
for ace in 1o*fasta.cap.ace; do
	perl ~/scripts/cap3/ace.ctg.reads.pl ${ace} > ${ace}.contig.kmers
done
