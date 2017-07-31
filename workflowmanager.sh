#!/bin/bash
#intended for linux/macOS

#USAGE: [-ref path_to_referencegenome] [-genomedir path_to_genome_directory] [-CPUS [num]]
#example
#mbair: bash workflowmanager.sh -ref ~/Downloads/Example1/Genomes/EEE_Florida91-4697.fasta -genomedir ~/Downloads/Example1/Genomes -CPUS 10

### SBATCH vars
#SBATCH --job-name=campylobacter-parsnpvsksnp-CPU-10-2GB
#
# Project:
#SBATCH --account=nn9305k
#
# Wall clock limit:
#SBATCH --time=02:20:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=2GB
#
# Number of cores:
#SBATCH --cpus-per-task=20
#SBATCH --output=slurmLogs/slurm%j.out
# Notifications
#SBATCH --mail-type=END,FAIL,BEGIN # notifications for job done & fail & started
#SBATCH --mail-user=sebasso@ifi.uio.no



###    log structure
# $currdir/"results"/$NOW-results"/stdout
# $currdir/"results"/$NOW-results"/stderr
# IF runned like this: bash workflowmanager.sh -ref ~/Genomes/EEE_Florida91-4697.fasta -genomedir ~/Genomes > stdout 2> stderr
# $currdir/"results"/$NOW-results"/$parsnp_output/parsnp.stdout
# $currdir/"results"/$NOW-results"/$parsnp_output/parsnp.stderr
# $currdir/"results"/$NOW-results"/$ksnp_output/ksnp.stdout
# $currdir/"results"/$NOW-results"/$ksnp_output/ksnp.stderr

source flows/cleanup.sh

trap cleanup 1 2 3 9 15

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
      -abel)
      shift
      abel="$1"
      shift
      ;;
      -name)
      shift
      jobname="$1"
      shift
      ;;
      -simulate)
      shift
      simulate="$1"
      shift
      ;;
      -percent)
      shift
      percent="$1"
      shift
      ;;
      -num_genomes)
      shift
      num_genomes="$1"
      shift
      ;;
      -probabilities)
      shift
      probabilities="$1"
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

if [ ! -z "$abel" ];
then
  source flows/abel.sh
else
  NOW=$(date +"%Y-%b-%d-%H:%M")
fi

#local variables setup
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


if [ ! -z "$simulate" ];
then
  if [ -z "$percent" ];
  then
    percent=2
  fi
  if [ -z "$num_genomes" ];
  then
    num_genomes=30
  fi
  printf "\ncurrdir: "$currdir"\n"
  simulated_dir=$currdir/evolve_genome/simulated_genomes/$NOW$jobname
  printf "simulated_dir: "$simulated_dir"\n"
  mkdir -p $simulated_dir
  genome_path=$simulated_dir/genomes
  mkdir $genome_path
  snp_stats=$simulated_dir/snp_stats
  mkdir $snp_stats

  if [ -z "$probabilities" ]; #TODO: RUN fasttree on $simulated_dir/fasttreeoutput.fasta and output it to $simulated_dir
  then
    python evolve_genome/genome_evolver.py --refgenome $ref --percent $percent --num_genomes $num_genomes --outputdir $simulated_dir
  else
    python evolve_genome/genome_evolver.py --refgenome $ref --percent $percent --num_genomes $num_genomes --outputdir $simulated_dir --nucleotide_probabilities $probabilities
  fi
fi

#### preprocessing of inputfiles(ksnp ONLY((Galaxy inputformat legacy))
#list files as one args as in galaxy framework when <param type=data multiple=true> , delimited set of inputfiles
input_files=""
for f in $genome_path/*
do
  input_files=$input_files","$f
done
input_files=${input_files:1} # remove first ,

ksnp_output="$kSNP_path/ksnp_output/$NOW"
parsnp_output="$parsnp_path/parsnp_output/$NOW"


# setup logfiles and output
mkdir -p $ksnp_output
mkdir -p $parsnp_output
touch $ksnp_output/stderr
touch $ksnp_output/stdout
touch $parsnp_output/stderr
touch $parsnp_output/stdout


# loading tools run defintions(source runs the code in this shell)
source flows/tool-config.sh
# SPAWNING two subshells, each for doing a different WGA- each fork costs 2ms ->  $! stores the PID of the LAST executed command
pids=() #for printing
tool_names=()
counter=0
for i in ${tool_array[@]};
do
   ( $i ) &  #calls tool_methods as a subshell -> enables tools parallelism
   pids[counter]=$!
   tool_names[counter]="$i"
   let counter=counter+1
done

printf "SHELL PID: $$\n"
printf '%s\t' "${tool_names[@]}"
printf '%s\t' "${pids[@]}"
printf "\n"

# waiting on tool subshells
counter=0
for pid in ${pids[@]}; #alt `jobs -p`but fails if processes finished before this
do
    printf "waiting on: ${tool_names[$counter]} $pid"
    wait $pid
    exit_status=$?
    if [[ $exit_status -ne 0 ]]; then
      >&2 printf "\ntool: ${tool_names[$counter]} $pid failed with status: $exit_status \n"
      cleanup
    fi
    printf "\tsuccess: ${tool_names[$counter]} $pid finished\n"
    let counter=counter+1
done
printf "Genome alignment done\n"
date



main_result_folder=$currdir/"results"
run_specific=$main_result_folder"/$NOW-results$jobname"

if [ ! -d  $main_result_folder ]; then
      mkdir $main_result_folder
fi

mkdir -p $run_specific
printf "\nrun specific: $run_specific\n"


if [ ! -z "$simulate" ];
then
  printf "Simulation started:\n"
  simulation_res=$run_specific/simulation # TODO: fix output
  #TODO: fix snp format from simulator
  printf "snp_statsdirectory:",$snp_stats
  ls $snp_stats
  printf "\n"
  mkdir -p $simulation_res/ksnp
  mkdir -p $simulation_res/parsnp
  python $currdir/snp_comparator.py $simulation_res/ksnp $ksnp_output/kSNP_SNPs_POS_formatted.tsv $snp_stats/reference_formatted_snps.tsv
  python $currdir/snp_comparator.py $simulation_res/parsnp $parsnp_output/parsnp_snps_sorted.tsv $snp_stats/reference_formatted_snps.tsv

  cd $currdir/common-sw
  if [ "$OS" == "Darwin" ];
  then
    printf "currdir: "$currdir
    ./FastTreeMP_osx -slow -nt -gtr $simulated_dir/fasttreeoutput.fasta > $simulated_dir/reference-tree.tree
  elif [ "$OS" == "Linux" ]; then
    ./FastTreeMP_linux -slow -nt -gtr $simulated_dir/fasttreeoutput.fasta > $simulated_dir/reference-tree.tree
  fi
  cd $currdir
  python $currdir/tree_comparator.py $simulation_res/parsnp $parsnp_output/parsnp.tree $simulated_dir/reference-tree.tree
  python $currdir/tree_comparator.py $simulation_res/ksnp $ksnp_output/ksnp.tree $simulated_dir/reference-tree.tree

  touch $simulation_res/mergedtrees.tree
  cat $ksnp_output/ksnp.tree >> $simulation_res/mergedtrees.tree
  printf "\n" >> $simulation_res/mergedtrees.tree
  cat $simulated_dir/reference-tree.tree >> $simulation_res/mergedtrees.tree
  java -jar $currdir/tree_comp/bin/TreeCmp.jar -m -d ms rf pd qt -i $simulation_res/mergedtrees.tree\
   -o $simulation_res/ksnp/distances-refvsparsnp.txt -I -P

  cat $parsnp_output/parsnp.tree > $simulation_res/mergedtrees.tree
  printf "\n" >> $simulation_res/mergedtrees.tree
  cat $simulated_dir/reference-tree.tree >> $simulation_res/mergedtrees.tree
  java -jar $currdir/tree_comp/bin/TreeCmp.jar -m -d ms rf pd qt -i $simulation_res/mergedtrees.tree\
   -o $simulation_res/parsnp/distances-refvsksnp.txt -I -P


  printf "after simulation res\n"
fi

##### COMPARATORS #######
##### SNP comparison
python $currdir/snp_comparator.py $run_specific $ksnp_output/kSNP_SNPs_POS_formatted.tsv $parsnp_output/parsnp_snps_sorted.tsv
printf "SNP comparison -> Done \n"
date
### ML tree comparison
#RF distance, common leaves, number of leaves, and total nodes.
python $currdir/tree_comparator.py $run_specific $ksnp_output/ksnp.tree $parsnp_output/parsnp.tree

# Matching Split distance, Robinson-Foulds distance, path difference distance, quartet distance
# number of taxan
touch $run_specific/mergedtrees.tree
cat $ksnp_output/ksnp.tree >> $run_specific/mergedtrees.tree
printf "\n" >> $run_specific/mergedtrees.tree
cat $parsnp_output/parsnp.tree >> $run_specific/mergedtrees.tree
java -jar $currdir/tree_comp/bin/TreeCmp.jar -m -d ms rf pd qt -i $run_specific/mergedtrees.tree\
 -o $run_specific/tree-distances.txt -I -P

printf "\nTree comparison -> Done \n"



###    Space saving    ### (moved here so there is less to cp below:)
#xmfa gets really large really fast, keep binary file parsnp.ggr | saves around 80 times disk space
# -> to get it back RUN: ./harvest_osx -i $run_specific/parsnp.ggr -X $run_specific/parsnp.xmfa #inside parsnp_folder

# Moving results to datefolder in results/
mv $ksnp_output/ksnp.tree $run_specific
mv $ksnp_output/kSNP_SNPs_POS_formatted.tsv $run_specific
mv $ksnp_output/SNPs_all $run_specific
mv $ksnp_output/SNPs_all_matrix.fasta $run_specific

mv $parsnp_output/parsnp.tree $run_specific
#mv $parsnp_output/parsnp_SNPs_POS_formatted.tsv $run_specific
mv $parsnp_output/parsnp.ggr $run_specific
mv $parsnp_output/parsnp.xmfa $run_specific
mv $parsnp_output/parsnp.snps.vcf $run_specific
mv $parsnp_output/parsnp_snps_sorted.tsv $run_specific

ksnp_logs=$run_specific/ksnp_logs
parsnp_logs=$run_specific/parsnp_logs
mkdir -p $ksnp_logs
mkdir -p $parsnp_logs

mv $ksnp_output/stderr $ksnp_logs
mv $ksnp_output/stdout $ksnp_logs

mv $parsnp_output/stderr $parsnp_logs
mv $parsnp_output/stdout $parsnp_logs

date
cleanup_junk
exit 0
