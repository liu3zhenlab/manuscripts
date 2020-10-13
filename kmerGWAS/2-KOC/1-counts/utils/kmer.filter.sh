### merge:
for txt in ../../7-kmers/result/*txt; do
	echo $txt
	geno=$(echo $txt | sed 's/.*\///g' | sed 's/.txt//g')
	out="SB."$geno.sbatch
	echo "#!/bin/bash -l" > $out
	echo "#SBATCH --mem-per-cpu=4G" >> $out
	echo "#SBATCH --time=0-23:00:00" >> $out
	echo "#SBATCH --nodes=1" >> $out
	echo "#SBATCH --ntasks-per-node=1" >> $out
	echo "#SBATCH --partition=ksu-plantpath-liu3zhen.q,batch.q,killable.q" >> $out
	echo "awk '\$2 >= 5' "$txt" > "$geno".txt" >> $out
	sbatch $out
done

