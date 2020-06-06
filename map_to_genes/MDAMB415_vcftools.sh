#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/genes
cd ./MDAMB415
conda activate vcftools
vcftools --gzvcf MDAMB415_vs_genes.vcf.gz --minQ 30 --minDP 3 --recode --recode-INFO-all --out MDAMB415_filtered.q30.mindp3.vcf
module load samtools
cp MDAMB415_filtered.q30.mindp3.vcf.recode.vcf MDAMB415_filtered.q30.mindp3.vcf.bk
bgzip MDAMB415_filtered.q30.mindp3.vcf.recode.vcf
tabix MDAMB415_filtered.q30.mindp3.vcf.recode.vcf.gz
/shelf/apps/pjt6/conda/envs/bcftools/bin/bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  MDAMB415_filtered.q30.mindp3.vcf.recode.vcf.gz > MDAMB415_H1_SNPs_reconstruct.q30.mindp3.fasta
