#!/bin/bash
#===============================================================================
# File Name    : s14_hmmsearch.sh
# Description  : Identify metabolic genes in isolate genomes
# Usage        : sbatch s04_hmmsearch.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2023-09-18
# Last Modified: 2023-09-18
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=HMMER
#SBATCH --array=1-501%35
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm/hmmsearch/x_hmmsearch_%A_%a.out
#SBATCH --error=slurm/hmmsearch/y_hmmsearch_%A_%a.out

#module load hmmer
eval $( spack load --sh hmmer )

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/bakta"
outdir="${basedir}/isolate_hmmsearch_KEGG"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/samplelist.txt`

set -x

time hmmsearch --tblout ${outdir}/${sample}_hmmsearch_KEGG.out -E 0.05 /ref/gdlab/software/envs/VIBRANT-Master/VIBRANT/databases/KEGG_profiles_prokaryotes.HMM ${indir}/${sample}/${sample}.faa

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