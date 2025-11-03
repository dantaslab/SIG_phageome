#!/bin/bash
#===============================================================================
#
# File Name    : s08b_Cenote-Taker2_sum.sh
# Description  : Summarize CT2 CONTIG_SUMMARY.tsv into a single file
# Usage        : s01b_Cenote-Taker2_sum.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 16 2021
# Modified On  : Jan 13 2023
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=sum
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/ct2_sum/x_ct2sum_%A.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/ct2_sum/y_ct2sum_%A.err

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/ct2"
samplelist="${basedir}/samplelist.txt"
outfile="${indir}/ct2_sum.tsv"

# create output file
touch ${outfile}

header="ORIGINAL_NAME\tCENOTE_NAME\tORGANISM_NAME\tEND_FEATURE\tLENGTH\tORF_CALLER\tNUM_HALLMARKS\tHALLMARK_NAMES\tBLASTN_INFO"
echo -e ${header} >> ${outfile}

# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        CT_Data=`grep -P "${sample}" ${indir}/${sample}/${sample}_CONTIG_SUMMARY.tsv`
        echo -e "${CT_Data}" >> ${outfile}

done < ${samplelist}
