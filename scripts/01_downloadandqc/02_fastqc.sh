#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

hostname
date

#################################################################
# FastQC
#################################################################
module load fastqc/0.12.1
module load parallel/20180122

# set input/output directory variables
INDIR=../../data/fastq
REPORTDIR=../../results/01_downloadandqc/fastqc_reports
mkdir -p $REPORTDIR

ACCLIST=../../metadata/SRR_Acc_List.txt

# run fastp in parallel
cat $ACCLIST | parallel -j 2 \
    "fastqc --outdir $REPORTDIR $INDIR/{}_{1..2}.fastq.gz"

date
