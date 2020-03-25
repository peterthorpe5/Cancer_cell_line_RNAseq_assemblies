#!/bin/bash -l
#SBATCH -J HCC_1395_RNAseq
#SBATCH --tasks-per-node=8
#SBATCH -p bigmem
#SBATCH --mem=100GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/
# make a folder of the cell name
mkdir HCC_1395
cd ./HCC_1395
#rm -rf HCC_1395
conda activate trinity2.9.1
#cp ../SRR925703_*.gz ./
conda activate trinity2.9.1

# QC trim the reads and remove vector
#trimmomatic PE -summary trim_summary.txt  -threads 8 -phred33 SRR925703_1.fastq.gz SRR925703_2.fastq.gz HCC_1395_paired_1.fq crap1.fastq.gz HCC_1395_paired_2.fq crap2.fastq.gz ILLUMINACLIP:~/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:47
#delete the raw reads
rm *.fastq.gz

# de nove RNAseq assembly: genome guided may not represent the cancer transcripts
#Trinity --full_cleanup --seqType fq --left HCC_1395_paired_1.fq --right HCC_1395_paired_2.fq --CPU 8 --max_memory 250G --min_glue 3 --output HCC_1395_trinity
#cp HCC_1395_trinity/Trinity.fasta HCC_1395.trinity.fasta

#QC the assembly and trash stuff .. this is a bit gung ho!!!
#transrate --assembly HCC_1395_trinity.Trinity.fasta --threads 8 --left HCC_1395_paired_1.fq --right HCC_1395_paired_2.fq

# transdecoder to extract the cds
#TransDecoder.LongOrfs -t HCC_1395_trinity.Trinity.fasta
#diamond blastp --query ./*/longest_orfs.pep --db /gpfs1/scratch/bioinf/db/databases/uniprot.dmnd --threads 8 --max-target-seqs 1 --outfmt 6 -o blastp.outfmt6
#TransDecoder.Predict -t HCC_1395_trinity.Trinity.fasta --retain_blastp_hits blastp.outfmt6

# GENOME guided assembly - map the reads
 #STAR --genomeDir ../star_indicies/ --limitGenomeGenerateRAM 5554136874 --limitBAMsortRAM 5554136874 --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix HCC_1395 --readFilesIn  HCC_1395_paired_1.fq HCC_1395_paired_2.fq
# GENOME guided assembly. Assemble from the mapped reads
 #Trinity --genome_guided_bam HCC_1395Aligned.sortedByCoord.out.bam --genome_guided_max_intron 10000 --max_memory 10G --CPU 12 --max_memory 100G --full_cleanup --output HCC_1395_GG_Trinity --genome_guided_min_coverage 5
cp ./HCC_1395_GG_Trinity/Trinity-GG.fasta ./HCC_1395_GG.fasta
TransDecoder.LongOrfs -t HCC_1395_GG.fasta
diamond blastp --query ./*/longest_orfs.pep --db /gpfs1/scratch/bioinf/db/databases/uniprot.dmnd --threads 8 --max-target-seqs 1 --outfmt 6 -o blastp.outfmt6
TransDecoder.Predict -t HCC_1395_GG.fasta --retain_blastp_hits blastp.outfmt6

