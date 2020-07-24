#!/bin/bash -l
#SBATCH -J HCC1937_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir HCC1937
cd ./HCC1937
conda activate seq2HLA
seq2HLA -1 ../HCC1937_paired_1.fq.gz -2 ../HCC1937_paired_2.fq.gz -r HCC1937_trimmed -p 16 
