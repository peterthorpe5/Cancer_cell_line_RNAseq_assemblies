from collections import defaultdict


info_table="""HCC_1395	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925703/SRR925703_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925703/SRR925703_2.fastq.gz
HCC_1806	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925708/SRR925708_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925708/SRR925708_2.fastq.gz
MDAMB361	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925727/SRR925727_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925727/SRR925727_2.fastq.gz
MDAMB453	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/008/SRR1283038/SRR1283038_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/008/SRR1283038/SRR1283038_2.fastq.gz
HCC1954	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925710/SRR925710_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925710/SRR925710_2.fastq.gz
MCF7	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925723/SRR925723_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925723/SRR925723_2.fastq.gz
MDAMB231	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925726/SRR925726_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925726/SRR925726_2.fastq.gz
MDAMB415	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR253/009/SRR2532369/SRR2532369_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR253/009/SRR2532369/SRR2532369_2.fastq.gz
SKBR3	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925729/SRR925729_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925729/SRR925729_2.fastq.gz
T47D	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925736/SRR925736_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925736/SRR925736_2.fastq.gz
ZR751	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925740/SRR925740_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925740/SRR925740_2.fastq.gz
ZR75B	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925742/SRR925742_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925742/SRR925742_2.fastq.gz
AU565	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925694/SRR925694_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925694/SRR925694_2.fastq.gz
HCC1937	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925709/SRR925709_1.fastq.gz	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR925/SRR925709/SRR925709_2.fastq.gz""".split("\n")


# note downloading file like this doesnt work in trinity. SO have to use SRAtool.
# Which wont work on Redhat6! FFS!!!
# S have to use a different cluster to donwload anc opy them over due to LIBC
# problems
# e.g.: /sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR925703

# dowload entries...
for entry in info_table:
    cell_line, R1, R2 = entry.split("\t")
    R1 = R1.rstrip()
    R2 = R2.rstrip()
    SRA = R1.split("/")[-1]
    SRA = SRA.split("_1")[0]
    cell_line = cell_line.rstrip()
    outsh = "%s_map_to_genes.sh" % cell_line.rstrip()
    # print("fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files %s\n" %SRA)

    f_out = open(outsh, "w")
    # 6 cores, so only 1 runs on a "node*" at a time. 
    print("\nqsub -q centos7.q -V -pe multi 2 %s_map_to_genes.sh\n" % cell_line)
    f_out.write("#!/bin/bash\n")
    f_out.write("#$ -cwd\n")
    f_out.write("#$ -pe multi 2\n")
    f_out.write("cd /storage/home/users/pjt6/cancer_cell_lines/map_to_genes\n")

    mkdir = "# make a folder of the cell name\nmkdir %s\n" % cell_line
    f_out.write(mkdir)
    cd_into = "cd ./%s\n" % cell_line
    f_out.write(cd_into)

    # q30 min len 57
    trimmo = "\n# QC trim the reads and remove vector\njava -jar /shelf/training/Trimmomatic-0.38/trimmomatic-0.38.jar PE -summary trim_summary.txt  -threads 4 -phred33 /storage/home/users/pjt6/cancer_cell_lines/%s_1.fastq.gz /storage/home/users/pjt6/cancer_cell_lines/%s_2.fastq.gz %s_paired_1.fq  crap1.fastq.gz %s_paired_2.fq crap2.fastq.gz  ILLUMINACLIP:/shelf/training/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:57\n" % (SRA, SRA, cell_line, cell_line)
    f_out.write(trimmo)
    delete_SRA_files = "#delete the raw reads\nrm crap*\n"
    f_out.write(delete_SRA_files)
    # trinity de novo cmd
    star_cmd = " /shelf/apps/pjt6/apps/STAR-master/bin/Linux_x86_64_static/STAR --genomeDir ../star_indicies/   --limitGenomeGenerateRAM 5554136874 --limitBAMsortRAM 5554136874 --runThreadN 4 --outSAMtype BAM SortedByCoordinate     --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix %s --readFilesIn  %s_paired_1.fq %s_paired_2.fq\n" % ( cell_line, cell_line, cell_line)
    f_out.write(star_cmd)
    index_s = "samtools index *.bam\n"
    f_out.write(index_s)
    # freebayes
    f_cmd = "freebayes -b %s*.bam -f ../Homo_sapiens.GRCh38.cds.all.fa  --no-population-priors --min-alternate-count 5 > %s_vs_genes.vcf\n" % (cell_line, cell_line)
    f_out.write(f_cmd)
    f_out.write("module load samtools\n")
    f_out.write("cp *.vcf vcf.bk\n")

    f_out.write("bgzip *.vcf\n")
    f_out.write("tabix *.vcf.gz\n")
    #now ready for bcftools to re make the genes

    # close this so we can open another file
    f_out.close()
    
