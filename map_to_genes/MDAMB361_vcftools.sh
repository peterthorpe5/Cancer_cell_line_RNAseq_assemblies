#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./MDAMB361
conda activate vcftools
vcftools --gzvcf MDAMB361_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out MDAMB361_filtered.q30.mindp3.vcf
module load samtools
cp MDAMB361_filtered.q30.mindp3.vcf.recode.vcf MDAMB361_filtered.q30.mindp3.vcf.bk
bgzip MDAMB361_filtered.q30.mindp3.vcf.recode.vcf
tabix MDAMB361_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  MDAMB361_filtered.q30.mindp3.vcf.recode.vcf.gz > MDAMB361_H1_SNPs_reconstruct.q30.mindp3.fasta
