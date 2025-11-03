#!/bin/bash
#===============================================================================
# File Name    : extract_AMG.sh
# Description  : Identify metabolic genes in isolate genomes
# Usage        : sbatch extract_AMG.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2023-09-18
# Last Modified: 2023-10-02
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=HMMER
#SBATCH --array=1-493%50
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm/hmmsearch/x_clean_%A_%a.out
#SBATCH --error=slurm/hmmsearch/y_clean_%A_%a.out


basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/isolate_hmmsearch_KEGG"
outdir="${basedir}/isolate_hmmsearch_KEGG_clean493"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/samplelist_493.txt`

set -x

sed '1,3d' ${indir}/${sample}_hmmsearch_KEGG.out | head -n -10 > ${outdir}/${sample}_hmmsearch_KEGG_clean.out

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