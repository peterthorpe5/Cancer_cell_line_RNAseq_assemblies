#!/bin/bash -l
#SBATCH -J ZR751_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir ZR751
cd ./ZR751
conda activate seq2HLA
seq2HLA -1 ../ZR751_paired_1.fq.gz -2 ../ZR751_paired_2.fq.gz -r ZR751_trimmed -p 16 
