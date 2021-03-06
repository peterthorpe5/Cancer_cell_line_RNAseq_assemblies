#!/bin/bash -l
#SBATCH -J HCC_1395_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir HCC_1395_Raw
cd ./HCC_1395_Raw
conda activate seq2HLA
seq2HLA -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925703_1.fastq.gz -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925703_2.fastq.gz -r HCC_1395_raw -p 16 
