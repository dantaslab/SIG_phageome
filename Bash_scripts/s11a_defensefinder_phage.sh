#!/bin/bash
#===============================================================================
# File Name    : s11_defensefinder.sh
# Description  : Annotation of assembled bacterial genomes, MAGs, and plasmids
# Usage        : sbatch s11_defensefinder.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Modified     : 2023-07-17
# Created      : 2023-07-19
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=defensefinder
#SBATCH --array=1-35
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/defensefinder/B_defensefinder_%a_%A.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/defensefinder/B_defensefinder_%a_%A.out

eval $( spack load --sh hmmer)

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/pharokka/change"
outdir="${basedir}/defensefinder_phage_rerun"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/phages/230721_DF_rerun.txt`

set -x
time defense-finder run ${indir}/${sample}/phanotate.faa -o ${outdir}/${sample} -w ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi