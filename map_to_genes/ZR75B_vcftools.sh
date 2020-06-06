#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./ZR75B
conda activate vcftools
vcftools --gzvcf ZR75B_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out ZR75B_filtered.q30.mindp3.vcf
module load samtools
cp ZR75B_filtered.q30.mindp3.vcf.recode.vcf ZR75B_filtered.q30.mindp3.vcf.bk
bgzip ZR75B_filtered.q30.mindp3.vcf.recode.vcf
tabix ZR75B_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  ZR75B_filtered.q30.mindp3.vcf.recode.vcf.gz > ZR75B_H1_SNPs_reconstruct.q30.mindp3.fasta
