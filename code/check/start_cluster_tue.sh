#!/bin/bash
#SBATCH --partition=mcs.default.q
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --mem=0
source $HOME/prepare_scip_dev.sh
./run_all.sh "$@"
