#!/bin/zsh

# USAGE: [runs(default 4)] [outputfile(default timestats)]
#example: zsh timemodule.sh 10 - runs the timing module ten time and writes output to timings/$NOW-resource-stats

printf "\n"
printf "@args:\n"
printf '%s\n' "$@"
#while cmd line args > 0 # ${#} = number of cmd arg lines -> https://learnxinyminutes.com/docs/bash/
while test ${#} -gt 0
do
  case "$1" in
      -ref)
      shift
      ref="$1"
      shift
      ;;
      -genomedir)
      shift
      genome_path="$1"
      shift
      ;;
      -CPUS)
      shift
      CPUS="$1"
      shift
      ;;
			-runs)
			shift
			runs="$1"
			shift
			;;
      *)
      echo "unknown cmd: $1"
      exit_module
      break
      ;;
  esac
done

if [ -z "$ref" ] && [ ! -f "$ref" ];
then
  printf "Reference genome must be supplied\n -ref path_to_referencegenome.fasta"
  exit 1
fi
if [ -z "$genome_path"  ] && [ ! -d "$genome_path" ];
then
  printf "Genome directory genome must be supplied\n -genome_dir path_to_genome_directory"
  exit 1
fi
if [ -z "$CPUS" ];
then
  CPUS=4
fi
if [ -z "$runs" ];
then
  runs=4
fi

#bash workflowmanager.sh \
#-ref /Users/sebastiansoberg/Downloads/Example1/Genomes/EEE_Florida91-4697.fasta  \
#-genomedir /Users/sebastiansoberg/Downloads/Example1/Genomes -CPUS 4 2>


OS=`uname`
printf "\nOS: $OS Cores: $CPUS\n"

opts='"%E\t%P\t%e"'
if [ $OS = "Darwin" ];
then
  timer=`which gtime`
	if [ -z $timer ]; then
		printf "Mac doesnt have gnu-time support, cannot format output\n"
    timer=`whereis time` # default timer without formatting option prints rusage
		exit 1
	fi
elif [ "$OS" = "Linux" ];
then
  timer=`which time`
  if [ -z $timer ]; then
      printf "unknown timer exiting."
      exit 1
  else
    timer=/usr/bin/time
  fi
else
	printf "\n OS not supported: $OS"
	exit 1
fi

printf "runs: $runs\n"
NOW=$(date +"%Y-%b-%d-%H:%M")
f1="timings/$NOW-resource-stats"
touch $f1
date
for i in {1..$runs}
do
	printf "Run $i\n"
	touch timings/tmp$NOW
	$timer -f $opts bash workflowmanager.sh -ref $ref -genomedir \
  $genome_path -CPUS $CPUS 2> timings/tmp$NOW > /dev/null
	printf "\nexit status: $?\n"
  tail -n 1 timings/tmp$NOW >> $f1
	printf "$tmp" >> $f1
	rm timings/tmp$NOW
	date
done


#eval $timer$opts bash workflowmanager.sh -ref ~/Downloads/Example1/Genomes/EEE_Florida91-4697.fasta -genomedir ~/Downloads/Example1/Genomes
printf "\nexit status: $?\n"
average=`cat $f1 | awk '{ total += $3 } END { print total/NR }'`
printf "Average time elapsed: \n" >> $f1
printf $average"\n" >> $f1
printf "\n*********************************************************\
\nJOB STATS available in:\n-->  $f1\
\n*********************************************************\n"
