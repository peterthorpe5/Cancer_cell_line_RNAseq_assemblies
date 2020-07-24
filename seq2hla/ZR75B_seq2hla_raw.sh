#!/bin/bash -l
#SBATCH -J ZR75B_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir ZR75B
cd ./ZR75B
conda activate seq2HLA
seq2HLA -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925742_1.fastq.gz -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925742_2.fastq.gz -r ZR75B_raw -p 16 
