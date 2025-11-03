#!/bin/bash
#===============================================================================
# File Name    : s04_vibrant.sh
# Description  : Virus Identification By iteRative ANnoTation
# Usage        : sbatch s05_vibrant.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 2.0
# Created On   : Nov 20 2021
# Last Modified: Jun 09 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=vibrant
#SBATCH --cpus-per-task=20
#SBATCH --mem=70G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/vibrant/x_vibrant_%A.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/vibrant/y_vibrant_%A.out

# In a Linux Screen
#srun -c 20 --time=20:00:00 --mem=70G --pty bash -i 
eval $( spack load --sh python@3.9.12 )
eval $( spack load --sh hmmer )
eval $( spack load --sh prodigal )
eval $( spack load --sh /knknaya ) # loading py-seaborn 
eval $( spack load --sh /yl6lkpn ) # loading py-matplotlib

python3 /ref/gdlab/software/envs/VIBRANT-Master/VIBRANT/VIBRANT_run.py -i /scratch/gdlab/nelson.m/SIG/phages/clean.fna -t 20 -virome -folder /scratch/gdlab/nelson.m/SIG/vibrant