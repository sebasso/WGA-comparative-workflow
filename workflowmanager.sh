#!/bin/bash
#intended use on linux server cluster


#USAGE: [-ref path_to_referencegenome] [-genomedir path_to_genome_directory] [-CPUS [num]]
#example bash workflowmanager.sh -ref ~/Downloads/Example1/Genomes/EEE_Florida91-4697.fasta -genomedir ~/Downloads/Example1/Genomes

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


exit_module(){
  printf "\nExiting"
  printf "["
  for i in {1..100}
  do
    sleep 0.002
    printf "#"
  done
  printf "]\n"
  printf "done\n"
  exit
}

printf 'cmd line args:\n%s\n' "$@"
while test ${#} -gt 0
do
  case "$1" in
      -ref)
      shift
      ref="$1"
      printf "\nref: $ref\n"
      printf "$1\n"
      shift
      ;;
      -genomedir)
      shift
      genome_path="$1"
      printf "\n gpath: $genome_path\n"
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

#optional option
if [ -z "$CPUS" ];
then
  CPUS=4
fi

#local variables
#NOW=$(date +"%Y-%b-%d-%H:%M")
currdir=`pwd`
kSNP_path=$currdir/ksnp
parsnp_path=$currdir

OS=`uname`
printf "OS: $OS\n "

#Setting up platform dependent executables for parsnp, ksnp/ is modified to fit any nix* platform
if [ "$OS" == "Darwin" ];
then
  parsnp_path="$parsnp_path/parsnp_mac"
elif [ "$OS" == "Linux" ];
then
  printf "$OS\n"
  parsnp_path="$parsnp_path/parsnp_linux"
else
  printf "\n Unknown os $OS"
  exit_module
fi

input_files=""
for f in $genome_path/*
do
  input_files=$input_files","$f
done
input_files=${input_files:1}

ksnp_output="$kSNP_path/ksnp_output"
parsnp_output="$parsnp_path/parsnp_output"

# spawning two subprocess, each for doing a WGS
cd $kSNP_path
( ./kSNP3 -in $input_files -outdir $ksnp_output -k 13 -CPU $CPUS -kchooser "1" -ML -path $kSNP_path ) &
cd ..
( $parsnp_path/./parsnp -r $ref -d $genome_path -o $parsnp_output -p $CPUS ) &

wait # wait for further execution on subprocesses spawned above
python $parsnp_path/parsnp_SNP_POS_extracter.py $parsnp_output

printf "\nGenome alignment are done\n"
python snp_comparator.py $kSNP_path/result_folder/kSNP_SNPs_POS_formatted.tsv $parsnp_output/parsnp_SNPs_POS_formatted.tsv
printf "\nComparison of genome alignment -> Done \n"
exit_module
