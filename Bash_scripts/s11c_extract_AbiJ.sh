#!/bin/bash
#===============================================================================
# File Name    : s11b_extract_sequnce.sh
# Description  : Extract gene sequences
# Usage        : sbatch s11b_extract_sequnce.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Modified     : 2023-10-12
# Created      : 2023-10-12
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=defensefinder
#SBATCH --array=1-105%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm/extract/x_AbiJ_%A_%a.out
#SBATCH --error=slurm/extract/y_AbiJ_%A_%a.out

basedir="/scratch/gdlab/nelson.m/SIG"
indir_defense="${basedir}/defensefinder"
indir_protein="${basedir}/bakta"
outdir="${basedir}/defense_gene/AbiJ"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/sp_AbiJ_list_105.txt`

set -x

proteinID=`grep 'AbiJ' ${indir_defense}/${sample}/${sample}_system_defense.tsv | cut -f7`
echo ">${sample}_${proteinID}" > ${outdir}/${sample}_AbiJ.fasta
grep "${proteinID}" -A1 ${indir_protein}/${sample}/${sample}.faa | grep -v ">${proteinID}" >> ${outdir}/${sample}_AbiJ.fasta

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