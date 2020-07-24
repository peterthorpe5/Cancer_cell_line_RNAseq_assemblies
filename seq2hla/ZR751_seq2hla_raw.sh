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
seq2HLA -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925740_1.fastq.gz -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925740_2.fastq.gz -r ZR751_raw -p 16 
