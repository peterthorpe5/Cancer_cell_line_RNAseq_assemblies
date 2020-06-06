#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./AU565
conda activate vcftools
vcftools --gzvcf AU565_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out AU565_filtered.q30.mindp3.vcf
module load samtools
cp AU565_filtered.q30.mindp3.vcf.recode.vcf AU565_filtered.q30.mindp3.vcf.bk
bgzip AU565_filtered.q30.mindp3.vcf.recode.vcf
tabix AU565_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  AU565_filtered.q30.mindp3.vcf.recode.vcf.gz > AU565_H1_SNPs_reconstruct.q30.mindp3.fasta
