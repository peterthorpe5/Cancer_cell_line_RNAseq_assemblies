#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/map_to_genes
cd ./HCC_1806
conda activate bcftools
bcftools consensus -s unknown -H 1 -f ../Homo_sapiens.GRCh38.cds.all.fa  HCC_1806_vs_genes.vcf.gz > HCC_1806_H1_SNPs_reconstruct.fasta
