#!/bin/bash
#SBATCH --job-name=cellranger2wk
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=60G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mail-user=user@uchc.edu
#SBATCH -o ../logs/%x_%j.out
#SBATCH -e ../logs/%x_%j.err

hostname
date

INDIR=../../data/fastq/
OUTDIR=../../results/02_cellrangercount/2wkold/
mkdir -p $OUTDIR

# rename fastqs to fit cellranger requirements
mv $INDIR/DRR608926_1.fastq.gz $INDIR/2wk_S1_R1_001.fastq.gz
mv $INDIR/DRR608926_2.fastq.gz $INDIR/2wk_S1_R2_001.fastq.gz

module load cellranger/8.0.1

cellranger count --id=2WKOLD \
--create-bam=false \
--include-introns=false \
--output-dir=$OUTDIR \
--transcriptome=../../../../SCRNA/mm10-2020-A_build/refdata-gex-mm10-2020-A/ \
--fastqs=$INDIR \
--sample=2wk \
--localcores=8 \
--localmem=60

date
