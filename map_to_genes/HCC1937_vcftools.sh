#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./HCC1937
conda activate vcftools
vcftools --gzvcf HCC1937_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out HCC1937_filtered.q30.mindp3.vcf
module load samtools
cp HCC1937_filtered.q30.mindp3.vcf.recode.vcf HCC1937_filtered.q30.mindp3.vcf.bk
bgzip HCC1937_filtered.q30.mindp3.vcf.recode.vcf
tabix HCC1937_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  HCC1937_filtered.q30.mindp3.vcf.recode.vcf.gz > HCC1937_H1_SNPs_reconstruct.q30.mindp3.fasta
