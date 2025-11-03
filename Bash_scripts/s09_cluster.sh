#!/bin/bash
#===============================================================================
#
# File Name    : s15_cluster.sh
# Description  : Cluster base on ANI
# Usage        : sbatch s08_cluster.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Nov 18 2021
# Modified On  : 2023-04-19
#===============================================================================
# Submission script for HTCF
#SBATCH --cpus-per-task=8
#SBATCH --job-name=blast
#SBATCH --mem=24G
#SBATCH --output=/scratch/gdlab/nelson.m/SIG/slurm/cluster/X_unbin_blast_95_nc_0.85_%a.out
#SBATCH --error=/scratch/gdlab/nelson.m/SIG/slurm/cluster/Y_unbin_blast_95_nc_0.85_%a.out


#module load ncbi-blast/2.6.0+
#eval $( spack load --sh blast-plus )
eval $( spack load --sh /u7ssbm4 )

basedir="/scratch/gdlab/nelson.m/SIG"
indir="${basedir}/phages"
outdir="${basedir}/sp_cluster"

# create a temp directory with outdir to deposit all the scaffold files
mkdir -p ${outdir}

makeblastdb -in ${indir}/sp_clean.fna -dbtype nucl -out ${outdir}/Cd_nt_db
blastn -query ${indir}/sp_clean.fna -db ${outdir}/Cd_nt_db -outfmt '6 std qlen slen' -max_target_seqs 100000 -out ${outdir}/unbin_blast_out.tsv -num_threads 8
python anicalc.py -i ${outdir}/unbin_blast_out.tsv -o ${outdir}/unbin_ani.tsv
python aniclust.py --fna ${indir}/sp_clean.fna --ani ${outdir}/unbin_ani.tsv --out ${outdir}/unbin_cluster.tsv --min_ani 95 --min_tcov 85 --min_qcov 0