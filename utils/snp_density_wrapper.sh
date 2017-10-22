#!/bin/bash


k=$1
directory=$2

score=0	
counter=0
printf "dir: $2\n"
printf "k: $1\n"
#printf "\n"
cd $2

files=`ls *.fna*`
#printf "files: $files"
for file in $files
do
	#printf "$file \n"
	num=$(python /Users/sebastiansoberg/WGA-comparative-workflow/utils/analyze_snp_density.py ${k} ${file})
	#printf "from analyzer: $num\n"
	counter=$(($counter+1))
	score=$(($score+num))
done

echo $score
echo $counter

result=$(($score/$counter))
printf "snps under k: $score\n avg: $result\n"

#cd -
