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
#SBATCH --array=1-501%50
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm/hmmsearch/x_KO_%A_%a.out
#SBATCH --error=slurm/hmmsearch/y_KO_%A_%a.out

#module load hmmer
eval $( spack load --sh hmmer )

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/isolate_hmmsearch_KEGG"
outdir="${basedir}/isolate_hmmsearch_KEGG"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/samplelist.txt`

set -x

grep 'K00248\|K00382\|K00390\|K00459\|K00558\|K00657\|K01047\|K01069\|K01220\|K01495\|K01737\|K01776\|K01915\|K01951\|K02788\|K03476\|K03800\|K06167\|K06920\|K10026\|K11212\|K15634\|K20895\|K22477' ${indir}/${sample}_hmmsearch_KEGG.out > ${outdir}/${sample}_AMG2.txt
awk 'BEGIN{OFS="\t"} {print $0, (FNR>1 ? FILENAME : "name")}' ${outdir}/${sample}_AMG2.txt > ${outdir}/${sample}_AMG.txt

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