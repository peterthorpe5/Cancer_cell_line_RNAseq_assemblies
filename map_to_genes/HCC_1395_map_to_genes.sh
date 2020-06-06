#!/bin/bash
#$ -cwd
#$ -pe multi 2
cd /storage/home/users/pjt6/cancer_cell_lines/map_to_genes
# make a folder of the cell name
mkdir HCC_1395
cd ./HCC_1395

# QC trim the reads and remove vector
java -jar /shelf/training/Trimmomatic-0.38/trimmomatic-0.38.jar PE -summary trim_summary.txt  -threads 4 -phred33 /storage/home/users/pjt6/cancer_cell_lines/SRR925703_1.fastq.gz /storage/home/users/pjt6/cancer_cell_lines/SRR925703_2.fastq.gz HCC_1395_paired_1.fq  crap1.fastq.gz HCC_1395_paired_2.fq crap2.fastq.gz  ILLUMINACLIP:/shelf/training/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:57
#delete the raw reads
rm crap*
 /shelf/apps/pjt6/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir ../star_indicies/   --limitGenomeGenerateRAM 5554136874 --limitBAMsortRAM 5554136874 --runThreadN 4 --outSAMtype BAM SortedByCoordinate     --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix HCC_1395 --readFilesIn  HCC_1395_paired_1.fq HCC_1395_paired_2.fq
samtools index *.bam
freebayes -b HCC_1395*.bam -f ../Homo_sapiens.GRCh38.cds.all.fa  --no-population-priors --min-alternate-count 5 > HCC_1395_vs_genes.vcf
module load samtools
cp *.vcf vcf.bk
bgzip *.vcf
tabix *.vcf.gz
