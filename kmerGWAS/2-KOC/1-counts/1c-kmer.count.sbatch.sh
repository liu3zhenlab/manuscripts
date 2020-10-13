#!/bin/bash
dataDir=../2-datamerge
outDir=1-counts
perl utils/kmer.sbatch.pl \
	--mem 3G --time 0-12:00:00 \
	--indir $dataDir --outdir $outDir \
	--fq1feature .R1.pair.fq.gz --fq2feature .R2.pair.fq.gz \
	--kmersize 25 \
	--mincount 1 \
	--hashsize "500M" \
	--fa2txt_script "utils/fasta2txt.pl" \
	--jellyfish "/homes/liu3zhen/local/bin/jellyfish" \
	--checkscript \
	--threads 8

