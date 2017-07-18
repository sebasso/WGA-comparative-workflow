#!/bin/bash

source /cluster/bin/jobsetup
module purge
module load python2

NOW=$(date +"%Y-%b-%d-%H:%M:%S")
