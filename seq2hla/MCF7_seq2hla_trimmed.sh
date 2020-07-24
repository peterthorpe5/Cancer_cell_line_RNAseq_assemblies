#!/bin/bash -l
#SBATCH -J MCF7_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir MCF7_trimmed
cd ./MCF7_trimmed
conda activate seq2HLA
seq2HLA -1 ../MCF7_paired_1.fq.gz -2 ../MCF7_paired_2.fq.gz -r MCF7_trimmed -p 16 
