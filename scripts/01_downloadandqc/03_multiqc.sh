#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

hostname
date

#################################################################
# Aggregate reports using MultiQC
#################################################################

module load MultiQC/1.15

INDIR=../../results/01_downloadandqc/fastqc_reports/
OUTDIR=../../results/01_downloadandqc/multiqc

# run on fastqc output
multiqc -f -o ${OUTDIR} ${INDIR}

date
