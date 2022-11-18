#!/bin/bash
# Job name:
#SBATCH --job-name=biotmle_update_deps
#
# Working directory:
#SBATCH --workdir=/global/home/users/nhejazi/
#
# Account:
#SBATCH --account=co_biostat
#
# Quality of Service:
#SBATCH --qos=biostat_savio2_normal
#
# Partition:
#SBATCH --partition=savio2
#
# Processors (1 node = 20 cores):
#SBATCH --nodes=1
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=12:00:00
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
export TMPDIR='~/rtmp'  # resolve update issues for compiled packages as per https://github.com/r-lib/devtools/issues/32
export R_LIBS_USER='/global/scratch/nhejazi/R'  # personal package library
module load gcc/6.3.0 r/3.6.3 r-packages/default
cd ~/biotmle-meta/simulations/
R CMD BATCH --no-save --no-restore \
  R/00_install_pkgs.R logs/00_install_pkgs.Rout
