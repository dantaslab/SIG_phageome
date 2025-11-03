#!/bin/bash

#===============================================================================
# File Name    : s10_panaroo_pharokka_pan.sh
# Description  : Runs the panaroo pangenome tool on genomes
# Usage        : sbatch s10_panaroo_pharokka_pan.sh
# Author       : Kailun Zhang
# Version      : 1.4
# Created On   : 2019-08-15 by Luke Diorio-Toth
# Last Modified: 2023-07-09
#===============================================================================

#SBATCH --job-name=panaroo
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/panaroo/x_panaroo_pan_%A.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/panaroo/y_panaroo_pan_%A.out

eval $( spack load --sh py-panaroo )

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/pharokka"
outdir="${basedir}/phage_panaroo"
tmpdir="${basedir}/tmp_${SLURM_JOB_ID}"

mkdir -p ${outdir}
mkdir -p ${tmpdir}

while read sample; do
    cp ${indir}/${sample}/pharokka.gff ${tmpdir}/${sample}.gff
done < ${basedir}/phages/finalphagelist.txt

set -x
time panaroo \
        -i ${tmpdir}/*.gff \
        -o ${outdir} \
        --clean-mode moderate \
        --core_threshold 0.1 \
        -c 0.5 \
        -a pan \
        --aligner mafft \
        -t ${SLURM_CPUS_PER_TASK}

RC=$?

rm -r ${tmpdir}

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi