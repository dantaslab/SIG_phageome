#!/bin/bash
#===============================================================================
# File Name    : s14_vcontact2.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s06_vcontact2.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.3
# Created On   : Dec 14 2021
# Last Modified: Apr 20 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=vcontact2
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/vcontact/x_vcontact2_%A.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/vcontact/y_vcontact2_%A.out


eval $( spack load --sh miniconda3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /home/kailun/conda/vContact2

#module load openjdk
eval $( spack load --sh /ubfljj4 )

basedir="/scratch/gdlab/nelson.m/SIG"

set -x
time vcontact2 --raw-proteins ${basedir}/checkv/tmp/proteins.faa --rel-mode 'Diamond' --proteins-fp ${basedir}/vcontact/g2g.csv --db 'ProkaryoticViralRefSeq207-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/kailun/conda/MAVERICLab-vcontact2-34ae9c466982/bin/cluster_one-1.0.jar --output-dir ${basedir}/vcontact/unbin


RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in vContact2"
  exit $RC
fi