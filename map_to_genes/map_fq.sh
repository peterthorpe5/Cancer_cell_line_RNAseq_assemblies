

trimmomatic PE -summary trim_summary.txt  -threads 8 -phred33 ../fq/SRR925703_1.fastq.gz ../fq/SRR925703_2.fastq.gz HCC_1395_paired_1.fq crap1.fastq.gz HCC_1395_paired_2.fq crap2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:47


 /shelf/apps/pjt6/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir ./star_indicies/ --limitGenomeGenerateRAM 5554136874 --limitBAMsortRAM 5554136874 --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix HCC_1395 --readFilesIn  HCC_1395_paired_1.fq HCC_1395_paired_2.fq

