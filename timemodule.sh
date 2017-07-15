#!/bin/zsh

# USAGE: [runs(default 4)] [outputfile(default timestats)]
#example: zsh timemodule.sh 10 - runs the timing module ten time and writes output to timings/$NOW-resource-stats

if [[ -z $1 ]]; then
	runs=4
else
	runs=$1
fi

printf "runs: $runs\n"
NOW=$(date +"%Y-%b-%d-%H:%M")
f1="timings/$NOW-resource-stats"
touch $f1
date
for i in {1..$runs}
do
	printf "Run $i\n"
	printf "$i: " >> $f1
	touch timings/tmp$NOW
	{ gtime -f "elapsed time: %E\tCPU: %P\tmax memory: %M (bytes)\n" bash workflowmanager.sh \
	-ref /Users/sebastiansoberg/Downloads/Example1/Genomes/EEE_Florida91-4697.fasta  \
	-genomedir /Users/sebastiansoberg/Downloads/Example1/Genomes -CPUS 4 2> timings/tmp$NOW > /dev/null; }
	tail -n 1 timings/tmp$NOW >> $f1
	printf "$tmp" >> $f1
	rm timings/tmp$NOW
	date
done

printf "\n*********************************************************\
\nJOB STATS available in:\n-->  $f1\
\n*********************************************************\n"

