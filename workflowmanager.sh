#!/bin/bash
#intended use on linux server cluster

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

printf '%s\n' "$@"
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


#NOW=$(date +"%Y-%b-%d-%H:%M")

#kSNP3/./kSNP3 -in in_list_scaffolds -outdir Res-campbyo/"$NOW-scaffolds" -k $kmer -CPU 10
currdir=`pwd`
kSNP_path=$currdir/ksnp
parsnp_path=$currdir


CPUS=10

OS=`uname`
if [ "$OS" == "Darwin" ];
then
  printf "$OS\n "
  #testing purposes
  parsnp_path="$parsnp_path/parsnp_mac"

  #TODO: change input format in KSNP main file
  #TODO: change input format parsnp
  #TODO: parsnp can put snp extracter outside main file
  input_files=""
  for f in $genome_path*
  do
    input_files="$input_files,$f"
  done

  ksnp_output="$kSNP_path/ksnp_output"
  parsnp_output="$parsnp_path/parsnp_output"
  printf "$kSNP_path\n"
  printf "$parsnp_path\n"
  file $kSNP_path/./kSNP3
  file $parsnp_path/./parsnp
  echo $CPUS
  echo $ref
  echo $genome_path
  echo $kSNP_path
  printf "ksnp output:\n$ksnp_output\n"
  printf "ksnp output:\n$parsnp_output\n"
  
  ( $kSNP_path/./kSNP3 -in $input_files -outdir $ksnp_output -k 13 -CPU $CPUS -kchooser "1" -ML -path $kSNP_path ) &
  #( $parsnp_path/./parsnp -r $ref -d $genome_path -o $parsnp_output -p $CPUS ) &

  wait
  printf "\ngenome alignment are done\n"
  #python snp_comparator.py $kSNP_path/$ksnp_output/kSNP_SNPs_POS_formatted.tsv $parsnp_path/$parsnp_output/kSNP_SNPs_POS_formatted.tsv
  printf "\nDone  with comparison of genome alignment\n"
  exit_module
elif [ "$OS" == "Linux" ];
then
  printf "$OS\n"
  parsnp_path="$parsnp_path/parsnp_linux"

  # ksnp -path is where ksnp executable is
  #change input_files on forehand to be on format filepath1,filepath2 from -d in parsnp
  input_files=""
  for f in $genome_path*
  do
    input_files="$input_files,$f"
  done

  ksnp_output="ksnp_output"
  parsnp_output="parsnp_output"
  ( $kSNP_path/.kSNP3 -in $input_files -outdir $ksnp_output -k 13 -CPU $CPUS -kchooser "1" -ML -path $kSNP_path ) &
  ( $parsnp_path./parsnp -r $ref -d $genome_path -o $parsnp_output -p $CPUS ) &

  wait
  printf "\ngenome alignment are done\n"
  python snp_comparator.py $kSNP_path/$ksnp_output/kSNP_SNPs_POS_formatted.tsv $parsnp_path/$parsnp_output/kSNP_SNPs_POS_formatted.tsv
  printf "\nDone  with comparison of genome alignment\n"
  exit_module
else
  printf "\n Unknown os $OS"
  exit_module
fi
