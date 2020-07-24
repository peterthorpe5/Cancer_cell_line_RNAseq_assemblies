#!/bin/bash -l
#SBATCH -J AU565_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir AU565
cd ./AU565
conda activate seq2HLA
seq2HLA -1 ../AU565_paired_1.fq.gz -2 ../AU565_paired_2.fq.gz -r AU565_trimmed -p 16 
