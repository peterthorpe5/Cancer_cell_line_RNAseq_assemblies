#!/bin/bash -l
#SBATCH -J T47D_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir T47D_trimmed
cd ./T47D_trimmed
conda activate seq2HLA
seq2HLA -1 ../T47D_paired_1.fq.gz -2 ../T47D_paired_2.fq.gz -r T47D_trimmed -p 16 
