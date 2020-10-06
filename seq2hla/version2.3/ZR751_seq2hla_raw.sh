#!/bin/bash -l
#SBATCH -J ZR751_RNAseq
#SBATCH --tasks-per-node=16
#SBATCH -p bigmem
#SBATCH --mem=60GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla
# make a folder of the cell name
mkdir ZR751_Raw
cd ./ZR751_Raw
conda activate seq2HLA
python /home/pjt6/scratch/cancer_genomes/seq2hla/seq2HLA/seq2HLA.py -1 /home/pjt6/scratch/cancer_genomes/fq/SRR925740_1.fastq -2 /home/pjt6/scratch/cancer_genomes/fq/SRR925740_2.fastq -r ZR751_raw -p 16 
