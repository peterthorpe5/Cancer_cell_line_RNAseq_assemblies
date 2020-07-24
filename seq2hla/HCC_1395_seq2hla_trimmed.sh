#!/bin/bash -l
#SBATCH -J HCC_1395_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir HCC_1395_trimmed
cd ./HCC_1395_trimmed
conda activate seq2HLA
seq2HLA -1 ../HCC_1395_paired_1.fq.gz -2 ../HCC_1395_paired_2.fq.gz -r HCC_1395_trimmed -p 16 
