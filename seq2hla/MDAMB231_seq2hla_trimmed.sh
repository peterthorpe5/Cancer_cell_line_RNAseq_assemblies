#!/bin/bash -l
#SBATCH -J MDAMB231_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir MDAMB231_trimmed
cd ./MDAMB231_trimmed
conda activate seq2HLA
seq2HLA -1 ../MDAMB231_paired_1.fq.gz -2 ../MDAMB231_paired_2.fq.gz -r MDAMB231_trimmed -p 16 
