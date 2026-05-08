#!/bin/bash
#SBATCH --job-name=downloadSRA
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

hostname
date

#################################################################
# Download fastq files from SRA
#################################################################

# load software
module load parallel/20180122
module load sratoolkit/3.0.1

# Data origin:
    # https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJDB18857

OUTDIR=../../data/fastq
    mkdir -p ${OUTDIR}


# 2 samples
ACCLIST=../../metadata/SRR_Acc_List.txt

# Set a custom temp directory
export TMPDIR=${OUTDIR}/tmp
mkdir -p ${TMPDIR}

# use parallel to download both accessions
cat $ACCLIST | parallel -j 2 --tmpdir ${TMPDIR} "fasterq-dump -O ${OUTDIR} {}"

# compress the files
parallel -j 2 --tmpdir ${TMPDIR} gzip ::: ${OUTDIR}/*.fastq

# Clean up temp directory
rm -rf ${TMPDIR}

date
