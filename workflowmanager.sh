#!/bin/bash
#intended for linux/macOS

#USAGE: [-ref path_to_referencegenome] [-genomedir path_to_genome_directory] [-CPUS [num]]
#example
#mbair: bash workflowmanager.sh -ref ~/Downloads/Example1/Genomes/EEE_Florida91-4697.fasta -genomedir ~/Downloads/Example1/Genomes


#SBATCH --job-name=parsnpvsksnp-CPU-10-2GB
#
# Project:
#SBATCH --account=nn9305k
#
# Wall clock limit:
#SBATCH --time=00:20:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=2GB
#
# Number of cores:
#SBATCH --cpus-per-task=10
#SBATCH --output=slurmLogs/slurm%j.txt


###    log structure
# $currdir/"results"/$NOW-results"/stdout
# $currdir/"results"/$NOW-results"/stderr
# IF runned like this: bash workflowmanager.sh -ref ~Genomes/EEE_Florida91-4697.fasta -genomedir ~/Genomes > stdout 2> stderr
# $currdir/"results"/$NOW-results"/$parsnp_output/parsnp.stdout
# $currdir/"results"/$NOW-results"/$parsnp_output/parsnp.stderr
# $currdir/"results"/$NOW-results"/$ksnp_output/ksnp.stdout
# $currdir/"results"/$NOW-results"/$ksnp_output/ksnp.stderr

exit_module(){
  printf "\nExiting"
  printf "["
  for i in {1..70}
  do
    sleep 0.003
    printf "#"
  done
  printf "]\n"
  exit
}


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
      *)
      echo "unknown cmd: $1"
      exit_module
      break
      ;;
  esac
done


#checking required options and type
if [ -z "$ref" ] && [ ! -f "$ref" ];
then
  printf "Reference genome must be supplied\n -ref path_to_referencegenome.fasta"
  exit_module
fi

if [ -z "$genome_path"  ] && [ ! -d "$genome_path" ];
then
  printf "Genome directory genome must be supplied\n -genome_dir path_to_genome_directory"
  exit_module
fi

#optional option - no paramater set default=4
if [ -z "$CPUS" ];
then
  CPUS=4
fi

#local variables setup
NOW=$(date +"%Y-%b-%d-%H:%M")
currdir=`pwd`
kSNP_path=$currdir/ksnp
parsnp_path=$currdir

OS=`uname`
printf "\nOS: $OS Cores: $CPUS\n"

#Setting up platform dependent executables for parsnp, ksnp/ is modified to fit any nix* platform
if [ "$OS" == "Darwin" ] || [ "$OS" == "Linux" ];
then
  parsnp_path="$parsnp_path/parsnp"
else
  printf "\n OS not supported: $OS"
  exit_module
fi

#### Galaxy inputformat: preprocessing of inputfiles
#list files as one args as in galaxy framework when <param type=data multiple=true> , delimited set of inputfiles
input_files=""
for f in $genome_path/*
do
  input_files=$input_files","$f
done
input_files=${input_files:1} # remove first ,

ksnp_output="$kSNP_path/ksnp_output"
parsnp_output="$parsnp_path/parsnp_output"


# setup logfiles
mkdir $ksnp_output
mkdir $parsnp_output
touch $ksnp_output/stderr
touch $ksnp_output/stdout
touch $parsnp_output/stderr
touch $parsnp_output/stdout

# loading tools run defintions
source tool-config.sh
# SPAWNING two subshells, each for doing a different WGA- each fork costs 2ms
## #   $! stores the PID of the LAST executed command
run_ksnp
ksnp_PID=$!

run_parsnp
parsnp_PID=$!
date

printf "SHELL PID: $$ parsnp pid: $parsnp_PID ksnp pid: $ksnp_PID\n"
printf "Waiting ..."


wait $parsnp_PID
printf "\n"
exit_status=$? #must be assigned or the variable will be lost after a simple if
stop=0
if [[ exit_status -ne 0 ]]; then
  stop=1
  >&2 printf "parsnp failed code: $exit_status \n"
fi
printf "parsnp done, waiting on ksnp\n"
date

wait $ksnp_PID
printf "\n"
exit_status=$?
if [[ exit_status -ne 0 ]]; then
    stop=1
    >&2 printf "ksnp failed code: $exit_status \n"
    exit_module
fi

if [[ $stop -ne 0 ]]; then
  >&2 printf "\n no comparison will be done due to the above problems\n"
  exit_module
fi
printf "ksnp done\n"
date



##### COMPARATORS #######
##### SNP comparison
printf "Genome alignment done\n"
python $currdir/snp_comparator.py $ksnp_output/kSNP_SNPs_POS_formatted.tsv $parsnp_output/parsnp_SNPs_POS_formatted.tsv
printf "SNP comparison -> Done \n"

### ML tree comparison
#not implemented -> when implemented run both as subshells with a wait





###    Space saving    ### (moved here so there is less to cp below:)
#xmfa gets really large really fast, keep binary file parsnp.ggr | saves around 80 times disk space
# -> to get it back RUN: ./harvest_osx -i $run_specific/parsnp.ggr -X $run_specific/parsnp.xmfa #inside parsnp_folder
rm $parsnp_output/parsnp.xmfa


# Moving results to datefolder in results/
main_result_folder=$currdir/"results"
run_specific=$main_result_folder"/$NOW-results"

if [ ! -d  $main_result_folder ]; then
      mkdir $main_result_folder
fi

mkdir $run_specific
mv snps_stats.json $run_specific
mv $ksnp_output $run_specific
mv $parsnp_output $run_specific

exit_module
