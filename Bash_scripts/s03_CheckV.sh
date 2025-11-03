#!/bin/bash
#===============================================================================
# File Name    : s03_CheckV.sh
# Description  : Assessing the quality of metagenome-assembled viral genomes
# Usage        : sbatch s03_CheckV.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 3.0
# Modified     : Nov 02 2021
# Created      : 2023-04-26
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkv
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/checkv/cx_cp_%A_%a.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/checkv/cy_cp_%A_%a.err

eval $( spack load --sh miniconda3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /ref/gdlab/software/envs/CheckV
export CHECKVDB=/ref/gdlab/software/envs/CheckV/checkv-db-v1.5


set -x
checkv end_to_end /scratch/gdlab/nelson.m/SIG/phages/clean.fna /scratch/gdlab/nelson.m/SIG/checkv/ -t 16

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
