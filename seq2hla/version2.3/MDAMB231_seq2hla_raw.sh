#!/bin/bash -l
#SBATCH -J MDAMB231_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir MDAMB231_Raw
cd ./MDAMB231_Raw
conda activate seq2HLA
python /home/pjt6/scratch/cancer_genomes/seq2hla/seq2HLA/seq2HLA.py -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925726_1.fastq -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925726_2.fastq -r MDAMB231_raw -p 16 
