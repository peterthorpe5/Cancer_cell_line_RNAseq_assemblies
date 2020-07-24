#!/bin/bash -l
#SBATCH -J MDAMB361_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir MDAMB361
cd ./MDAMB361
conda activate seq2HLA
seq2HLA -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925727_1.fastq.gz -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925727_2.fastq.gz -r MDAMB361_raw -p 16 
