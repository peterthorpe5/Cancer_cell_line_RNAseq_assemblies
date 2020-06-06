#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./ZR751
conda activate vcftools
vcftools --gzvcf ZR751_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out ZR751_filtered.q30.mindp3.vcf
module load samtools
cp ZR751_filtered.q30.mindp3.vcf.recode.vcf ZR751_filtered.q30.mindp3.vcf.bk
bgzip ZR751_filtered.q30.mindp3.vcf.recode.vcf
tabix ZR751_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  ZR751_filtered.q30.mindp3.vcf.recode.vcf.gz > ZR751_H1_SNPs_reconstruct.q30.mindp3.fasta
