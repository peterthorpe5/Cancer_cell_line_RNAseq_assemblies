#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./MDAMB453
conda activate vcftools
vcftools --gzvcf MDAMB453_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out MDAMB453_filtered.q30.mindp3.vcf
module load samtools
cp MDAMB453_filtered.q30.mindp3.vcf.recode.vcf MDAMB453_filtered.q30.mindp3.vcf.bk
bgzip MDAMB453_filtered.q30.mindp3.vcf.recode.vcf
tabix MDAMB453_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  MDAMB453_filtered.q30.mindp3.vcf.recode.vcf.gz > MDAMB453_H1_SNPs_reconstruct.q30.mindp3.fasta
