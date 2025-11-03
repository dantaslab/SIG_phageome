#!/bin/bash
#===============================================================================
# File Name    : s11b_extract_sequnce.sh
# Description  : Extract gene sequences
# Usage        : sbatch s11b_extract_sequnce.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Modified     : 2023-10-12
# Created      : 2023-10-16
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=defensefinder
#SBATCH --array=1-110%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm/extract/x_Gabija_%A_%a.out
#SBATCH --error=slurm/extract/y_Gabija_%A_%a.out

basedir="/scratch/gdlab/nelson.m/SIG"
indir_defense="${basedir}/defensefinder"
indir_protein="${basedir}/bakta"
outdir="${basedir}/defense_gene/Gabija"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/sp_Gabija_list_110.txt`
FJIINC_11240
set -x

proteinID1=`grep 'Gabija' ${indir_defense}/${sample}/${sample}_system_defense.tsv | cut -f7 | cut -d ',' -f1`
proteinID2=`grep 'Gabija' ${indir_defense}/${sample}/${sample}_system_defense.tsv | cut -f7 | cut -d ',' -f2`
profile1=`grep 'Gabija' ${indir_defense}/${sample}/${sample}_system_defense.tsv | cut -f9 | cut -d ',' -f1`
profile2=`grep 'Gabija' ${indir_defense}/${sample}/${sample}_system_defense.tsv | cut -f9 | cut -d ',' -f2`

echo ">${sample}_${proteinID1}_${profile1}" > ${outdir}/${sample}_Gabija.fasta
grep "${proteinID1}" -A1 ${indir_protein}/${sample}/${sample}.faa | grep -v ">" >> ${outdir}/${sample}_Gabija.fasta
echo ">${sample}_${proteinID2}_${profile2}" >> ${outdir}/${sample}_Gabija.fasta
grep "${proteinID2}" -A1 ${indir_protein}/${sample}/${sample}.faa | grep -v ">" >> ${outdir}/${sample}_Gabija.fasta

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