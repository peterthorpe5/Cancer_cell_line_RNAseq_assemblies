#!/bin/bash -l
#SBATCH -J HCC1954_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir HCC1954_trimmed
cd ./HCC1954_trimmed
conda activate seq2HLA
seq2HLA -1 ../HCC1954_paired_1.fq.gz -2 ../HCC1954_paired_2.fq.gz -r HCC1954_trimmed -p 16 
