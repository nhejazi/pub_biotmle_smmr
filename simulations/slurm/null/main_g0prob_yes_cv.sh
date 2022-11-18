#!/bin/bash
# Job name:
#SBATCH --job-name=cv_biotmle_sim_global_null_main_yes_g0prob
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
  '--args sim_type=global_null cv_folds=2 sl_spec=main positivity_issues=yes' \
  R/03_run_sim_sl.R logs/null/03d_g0prob_yes_sl_main_cv.Rout
