#! /bin/bash

Help() {
	echo "Arguments:"
	echo "[hisat_index_prefix] [reads_F.fasta.gz] [reads_R.fasta.gz]"
	echo
	echo "Options:"
	echo "h	print this help"
	echo "p	the prefix for the output alignment files"
	echo
	}

	

while getopts ":htp:" option
do
	
	case $option in
	h)
		Help
		exit;;
	p)
		p=$OPTARG
		shift;;
	esac
done

R1="$3"

R2="$4"

output="$p.bam"

echo "aligning $R1 and $R2 reads against $2 genome"
echo "Output file will be $output"


hisat2 -x $2 -1 $R1 -2 $R2 --max-intronlen 10000 --summary-file "${p}_hisat2_summary.txt" | samtools view -b | samtools sort -o "$output"

samtools index "$output"