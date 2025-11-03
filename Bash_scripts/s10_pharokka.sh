#!/bin/bash
#===============================================================================
# File Name    : s13_pharokka.sh
# Description  : This script will run prokka in parallel
# Usage        : sbatch s13_pharokka.sh
# Author       : Kailun Zhang. kailun@wustl.edu
# Version      : 1.1 (v1.3.0)
# Modified     : 2022-07-20 
# Created      : 2023-07-07
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=prokka
#SBATCH --array=1-1105%50
#SBATCH --mem=4G
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm/pharokka/x_pharokka_%A_%a.out
#SBATCH --error=slurm/pharokka/y_pharokka_%A_%a.err

#eval $( spack load --sh miniconda3 )
#CONDA_BASE=$(conda info --base)
#source $CONDA_BASE/etc/profile.d/conda.sh
#conda activate pharokka_env_1.3.0

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/phages/single"
outdir="${basedir}/pharokka"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/phages/finalphagelist.txt`

set -x
time python /ref/gdlab/software/envs/pharokka/v1.3.0/pharokka/bin/pharokka.py -i ${indir}/${sample}.fna -o ${outdir}/${sample} -d /ref/gdlab/software/envs/pharokka/v1.3.0/pharokka/databases -t ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi