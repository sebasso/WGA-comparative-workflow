#!/bin/bash
#$1 path to read files in fastq
#$2 path to outputdirectory for fasta files from fastq

for f in $1/*
do
	filename=`basename $f`
	printf "\nbasenamefile: $filename"
	filename="${filename%.*}" #remove .gz
	filename="${filename%.*}" #remove .fastq
	filename="$filename.fasta"
	printf "\n filename $filename \n"
	seqtk seq -a -q20 -n N $f > $2$filename
done
