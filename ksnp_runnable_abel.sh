#!/bin/bash

### SBATCH vars
#SBATCH --job-name=KSNP-campylobacter-CPU-6-10GB
#
# Project:
#SBATCH --account=nn9305k
#
# Wall clock limit:
#SBATCH --time=04:20:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=10GB
#
# Number of cores:
#SBATCH --cpus-per-task=6
#SBATCH --output=ksnpslurmLogs/slurm%j.out
# Notifications
#SBATCH --mail-type=END,FAIL,BEGIN # notifications for job done & fail & started
#SBATCH --mail-user=sebasso@ifi.uio.no

source /cluster/bin/jobsetup
module purge
module load python2


genomedir=/work/projects/nn9305k/compgenomes/rawdata/campy_reads
input_files=""
for f in $genomedir/*
do
  input_files=$input_files","$f
done
input_files=${input_files:1} # remove first ,

currdir=`pwd`
echo `pwd`
ksnp_output=$currdir/resultsKSNP
kSNP_path=$currdir/ksnp

#touch $ksnp_output/stdout
#touch $ksnp_output/stderr

cd $kSNP_path
./kSNP3 -in $input_files -outdir $ksnp_output -k 13 -CPU 12 \
-kchooser "1" -ML -path $kSNP_path  > $ksnp_output/stdout  2> $ksnp_output/stderr