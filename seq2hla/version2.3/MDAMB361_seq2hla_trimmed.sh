#!/bin/bash -l
#SBATCH -J MDAMB361_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir MDAMB361_trimmed
cd ./MDAMB361_trimmed
conda activate seq2HLA
python /home/pjt6/scratch/cancer_genomes/seq2hla/seq2HLA/seq2HLA.py -1 /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla/MDAMB361_paired_1.fq -2 /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla/MDAMB361_paired_2.fq -r MDAMB361_trimmed -p 16 
