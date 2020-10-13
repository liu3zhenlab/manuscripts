### merge:
for acgt in [ACGT]*; do
	echo $acgt
	l3=$(echo $acgt | sed 's/.*\///g')
	out="SB."$l3.sbatch
	echo "#!/bin/bash -l" > $out
	echo "#SBATCH --mem-per-cpu=124G" >> $out
	echo "#SBATCH --time=0-23:00:00" >> $out
	echo "cat "$acgt"/* | awk 'seen[\$1]++==10' > "$l3".merge" >> $out
	sbatch $out
done

