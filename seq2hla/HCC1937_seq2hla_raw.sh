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
seq2HLA -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925709_1.fastq.gz -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925709_2.fastq.gz -r HCC1937_raw -p 16 
