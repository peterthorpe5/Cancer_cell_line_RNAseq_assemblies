#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./MCF7
conda activate vcftools
vcftools --gzvcf MCF7_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out MCF7_filtered.q30.mindp3.vcf
module load samtools
cp MCF7_filtered.q30.mindp3.vcf.recode.vcf MCF7_filtered.q30.mindp3.vcf.bk
bgzip MCF7_filtered.q30.mindp3.vcf.recode.vcf
tabix MCF7_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  MCF7_filtered.q30.mindp3.vcf.recode.vcf.gz > MCF7_H1_SNPs_reconstruct.q30.mindp3.fasta
