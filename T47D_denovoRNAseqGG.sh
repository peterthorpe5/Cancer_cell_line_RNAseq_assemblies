#!/bin/bash -l
#SBATCH -J T47D_RNAseq
#SBATCH --tasks-per-node=8
#SBATCH -p bigmem
#SBATCH --mem=100GB
cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/
# make a folder of the cell name
mkdir T47D
cd ./T47D
#rm -rf T47D
conda activate trinity2.9.1
#cp ../SRR925736_*.gz ./
conda activate trinity2.9.1

# QC trim the reads and remove vector
#trimmomatic PE -summary trim_summary.txt  -threads 8 -phred33 SRR925736_1.fastq.gz SRR925736_2.fastq.gz T47D_paired_1.fq crap1.fastq.gz T47D_paired_2.fq crap2.fastq.gz ILLUMINACLIP:~/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:47
#delete the raw reads
rm *.fastq.gz

# de nove RNAseq assembly: genome guided may not represent the cancer transcripts
#Trinity --full_cleanup --seqType fq --left T47D_paired_1.fq --right T47D_paired_2.fq --CPU 8 --max_memory 250G --min_glue 3 --output T47D_trinity
#cp T47D_trinity/Trinity.fasta T47D.trinity.fasta

#QC the assembly and trash stuff .. this is a bit gung ho!!!
#transrate --assembly T47D_trinity.Trinity.fasta --threads 8 --left T47D_paired_1.fq --right T47D_paired_2.fq

# transdecoder to extract the cds
#TransDecoder.LongOrfs -t T47D_trinity.Trinity.fasta
#diamond blastp --query ./*/longest_orfs.pep --db /gpfs1/scratch/bioinf/db/databases/uniprot.dmnd --threads 8 --max-target-seqs 1 --outfmt 6 -o blastp.outfmt6
#TransDecoder.Predict -t T47D_trinity.Trinity.fasta --retain_blastp_hits blastp.outfmt6

# GENOME guided assembly - map the reads
 #STAR --genomeDir ../star_indicies/ --limitGenomeGenerateRAM 5554136874 --limitBAMsortRAM 5554136874 --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix T47D --readFilesIn  T47D_paired_1.fq T47D_paired_2.fq
# GENOME guided assembly. Assemble from the mapped reads
 #Trinity --genome_guided_bam T47DAligned.sortedByCoord.out.bam --genome_guided_max_intron 10000 --max_memory 10G --CPU 12 --max_memory 100G --full_cleanup --output T47D_GG_Trinity --genome_guided_min_coverage 5
cp ./T47D_GG_Trinity/Trinity-GG.fasta ./T47D_GG.fasta
TransDecoder.LongOrfs -t T47D_GG.fasta
diamond blastp --query ./*/longest_orfs.pep --db /gpfs1/scratch/bioinf/db/databases/uniprot.dmnd --threads 8 --max-target-seqs 1 --outfmt 6 -o blastp.outfmt6
TransDecoder.Predict -t T47D_GG.fasta --retain_blastp_hits blastp.outfmt6

