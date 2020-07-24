#!/bin/bash -l
#SBATCH -J T47D_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir T47D
cd ./T47D
conda activate seq2HLA
seq2HLA -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925736_1.fastq.gz -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925736_2.fastq.gz -r T47D_raw -p 16 
