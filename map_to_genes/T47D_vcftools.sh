#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./T47D
conda activate vcftools
vcftools --gzvcf T47D_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out T47D_filtered.q30.mindp3.vcf
module load samtools
cp T47D_filtered.q30.mindp3.vcf.recode.vcf T47D_filtered.q30.mindp3.vcf.bk
bgzip T47D_filtered.q30.mindp3.vcf.recode.vcf
tabix T47D_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  T47D_filtered.q30.mindp3.vcf.recode.vcf.gz > T47D_H1_SNPs_reconstruct.q30.mindp3.fasta
