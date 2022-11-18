#!/bin/bash
# Job name:
#SBATCH --job-name=biotmle_sim_mod_findings_main_small_effect
#
# Working directory:
#SBATCH --workdir=/global/home/users/nhejazi/
#
# Account:
#SBATCH --account=co_biostat
#
# Partition:
#SBATCH --partition=savio2
#
# Quality of Service:
#SBATCH --qos=biostat_savio2_normal
#
# Processors (1 node = 20 cores):
#SBATCH --nodes=1
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=336:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=nhejazi@berkeley.edu
#
# Job output:
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out
#
## Command(s) to run:
export TMPDIR='/global/scratch/nhejazi/rtmp'  # resolve update issues for compiled packages as per https://github.com/r-lib/devtools/issues/32
export R_LIBS_USER='/global/scratch/nhejazi/R'  # personal package library
module load gcc/6.3.0 r/3.6.3 r-packages/default
cd ~/biotmle-meta/simulations/

R CMD BATCH --no-save --no-restore \
  '--args sim_type=true_findings cv_folds=1 sl_spec=main prop_genes=0.3 effect_size=small' \
  R/03_run_sim_sl.R logs/pmod/03e_effect_small_sl_main.Rout
