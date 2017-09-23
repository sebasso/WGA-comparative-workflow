#!/bin/bash
#$1 path to read files in fastq
#$2 path to outputdirectory for fasta files from fastq

for f in $1/*
do
	filename=`basename $f`
	filename="${filename%.*}" #remove .fastq
	filename="$filename.fasta"
	printf "\n filename $filename"
	echo $2$filename
	seqtk seq -a -q20 -n N $f > $2$filename
done
