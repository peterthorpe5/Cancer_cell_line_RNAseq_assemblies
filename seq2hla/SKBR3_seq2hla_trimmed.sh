#!/bin/bash -l
#SBATCH -J SKBR3_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir SKBR3_trimmed
cd ./SKBR3_trimmed
conda activate seq2HLA
seq2HLA -1 ../SKBR3_paired_1.fq.gz -2 ../SKBR3_paired_2.fq.gz -r SKBR3_trimmed -p 16 
