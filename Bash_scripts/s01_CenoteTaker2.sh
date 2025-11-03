#!/bin/bash
#===============================================================================
# File Name    : s01_CenoteTaker2.sh
# Description  : Discover and annoate virus sequences/genomes
# Usage        : sbatch s01_CenoteTaker2.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 3.0
# Created      : Jan 09 2023 based on Jian and Luke
# Modified     : Jan 09 2023
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=cenote-taker2
#SBATCH --cpus-per-task=2
#SBATCH --array=1-501%15
#SBATCH --mem=32G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/ct2/x_ct2_%A_%a.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/ct2/y_ct2_%A_%a.err

# Use 'pyfasta split -n 10 input.fasta' to split the .fasta file if it's too big and the result files are suggested to be around 2-10Mbp
# Run 1-2 samples first and use 'seff <job ID>' to check the usage of cpus and mem. Then, adjust these two parameters for the rest jobs
# Generate a directory and put this script in it before running. 

eval $( spack load --sh miniconda3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate cenote-taker2_env
CENOTE_BASE="/ref/gdlab/software/envs/Cenote-Taker2"

# Set directories of the work
DATA_BASE="/scratch/gdlab/nelson.m/SIG"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${DATA_BASE}/samplelist.txt`
#sample_out=`sed -n ${SLURM_ARRAY_TASK_ID}p ${DATA_BASE}/230109_split_list_new.txt`
#mkdir /tmp/cenote

set -x

# cp ${DATA_BASE}/d03.2_unicycler/${sample}/assembly.fasta ${DATA_BASE}/d03.2_unicycler/${sample}/${sample}_assembly.fasta
# -r flag: Name of this run. A directory of this name will be created. Must be unique from older runs or older run will be renamed. Must be less than 18 characters, using ONLY letters, numbers and underscores (_)
# Thus I change all special charactors to underscore in the out files as 'sample_out'
time python ${CENOTE_BASE}/run_cenote-taker2.py \
    -c ${DATA_BASE}/assemblies/${sample}.fna \
    -r ${sample} \
    -p True \
    -m 32 \
    -t ${SLURM_CPUS_PER_TASK}
    #--scratch_directory /tmp/kailun_cenote
#    --known_strains blast_knowns \
#    --blastn_db /scratch/ref/gdlab/blast_db/nt_2022_03_08/nt \
    
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi