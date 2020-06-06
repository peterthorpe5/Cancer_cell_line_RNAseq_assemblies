#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./SKBR3
conda activate vcftools
vcftools --gzvcf SKBR3_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out SKBR3_filtered.q30.mindp3.vcf
module load samtools
cp SKBR3_filtered.q30.mindp3.vcf.recode.vcf SKBR3_filtered.q30.mindp3.vcf.bk
bgzip SKBR3_filtered.q30.mindp3.vcf.recode.vcf
tabix SKBR3_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  SKBR3_filtered.q30.mindp3.vcf.recode.vcf.gz > SKBR3_H1_SNPs_reconstruct.q30.mindp3.fasta
