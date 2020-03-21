db=<Reference & GTF path>
fqpath=<reads path>
pair1feature=<reads_1 suffix>
pair2feature=<reads_2 suffix>
threads=20
for fq1 in $fqpath/*$pair1feature; do ## search for files with .R1.fq.gz in fqpath as fq1 input
        echo $fq1 ## print fq1
        fq2=$(echo $fq1 | sed "s/$pair1feature/$pair2feature/g") ## replace .R1.fq.gz to .R2.fq.gz to form fq2 input
        echo $fq2 ## print fq2
        out=$(echo $fq1 | sed 's/.*\///g' | sed "s/$pair1feature//g") 
        ## sed 's/[.*][\/]/[]/g' replaces words start with anyone character and filled with any other characters until to a '/' to null 
        ## sed "s/$pair1feature//g" replaces '.R1.fq.gz' to null
        echo $out
        STAR --runThreadN $threads \
                --genomeDir $db \
                --readFilesIn $fq1 $fq2 \
                --readFilesCommand zcat \ ## gzip reads input
                --alignIntronMax 100000 \
                --alignMatesGapMax 100000 \
                --outFileNamePrefix $out \
                --outSAMattrIHstart 0 \
                --outSAMmultNmax 1 \
                --outSAMstrandField intronMotif \
                --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode GeneCounts \
                --outFilterMismatchNmax 2 \ ## mismatch 2 for B73 lines
                --outFilterMismatchNoverLmax 0.02 \ ## 2 mismatch/100 read length = 0.02
                --outFilterMatchNmin 50 \
                --outSJfilterReads Unique \
                --outFilterMultimapNmax 1 \
                --outSAMmapqUnique 60 \
                --outFilterMultimapScoreRange 2
done
