#!/bin/bash -l
# Job name:
#SBATCH --job-name=$JOBNAME
#
# Partition:
#SBATCH --partition=vulcan
#
# Account:
#SBATCH --account=vulcan
#
# Wall clock limit:
#SBATCH --time=2:00:00
#
# Processors
#SBATCH --ntasks=8
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=$EMAIL
#
## Run command

cd $SLURM_SUBMIT_DIR;
module load python
source activate conda36
python DoGenfire.py $GENFIRESTRING 
