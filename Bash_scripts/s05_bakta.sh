#!/bin/bash
#===============================================================================
# File Name    : bakta.sh
# Description  : Annotation of assembled bacterial genomes, MAGs, and plasmids
# Usage        : sbatch bakta.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Modified     : Fri Oct 14 08:14:59 CDT 2022
# Created      : Fri Oct 14 08:14:59 CDT 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=bakta
#SBATCH --array=1-501%15
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/bakta/x_bakta_%a_%A.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/bakta/c_bakta_%a_%A.out

eval $( spack load --sh py-bakta )

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/assemblies"
outdir="${basedir}/bakta"
dbdir="/ref/gdlab/data/bakta_db/db"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/samplelist.txt`

set -x
time bakta \
  --db ${dbdir} \
  --min-contig-length 200 \
  --prefix ${sample} \
  --output ${outdir}/${sample}.fasta \
  --threads ${SLURM_CPUS_PER_TASK} \
  --output ${outdir}/${sample} \
  ${indir}/${sample}.fna
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