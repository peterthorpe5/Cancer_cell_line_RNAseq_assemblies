#!/bin/bash -l
#SBATCH -J MDAMB415_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir MDAMB415_trimmed
cd ./MDAMB415_trimmed
conda activate seq2HLA
seq2HLA -1 ../MDAMB415_paired_1.fq.gz -2 ../MDAMB415_paired_2.fq.gz -r MDAMB415_trimmed -p 16 
