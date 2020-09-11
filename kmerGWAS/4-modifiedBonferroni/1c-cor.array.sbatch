#!/bin/bash
#SBATCH --array=1-1000
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
module load R
arrayid=$SLURM_ARRAY_TASK_ID
declare -i first=($arrayid-1)*1000+1
declare -i end=$arrayid*1000
echo $first $end
Rscript 1c-num.cor.gte.9.R $first $end

