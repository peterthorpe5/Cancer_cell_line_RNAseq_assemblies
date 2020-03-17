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


dowwnload = "download.sh"
dwnload = open(dowwnload, "w")

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
    outsh = "%s_denovoRNAseqGG.sh" % cell_line.rstrip()
    #print("fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files %s\n" %SRA)

    f_out = open(outsh, "w")
    # 6 cores, so only 1 runs on a "node*" at a time. 
    print("\nsbatch %s_denovoRNAseqGG.sh\n" % cell_line)
    f_out.write("#!/bin/bash -l\n")
    f_out.write("#SBATCH -J %s_RNAseq\n" % cell_line)
    f_out.write("#SBATCH --tasks-per-node=8\n")
    f_out.write("#SBATCH -p bigmem\n")
    f_out.write("#SBATCH --mem=100GB\n")
    f_out.write("cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/\n")

    mkdir = "# make a folder of the cell name\nmkdir %s\n" % cell_line
    f_out.write(mkdir)
    cd_into = "cd ./%s\n" % cell_line
    f_out.write(cd_into)
    rm_old_tri = "#rm -rf %s\n" % cell_line
    f_out.write(rm_old_tri)
    cmd1 = "\n# download the raw reads\nwget %s\n" % R1.rstrip()
    cmd2 = "wget %s\n" % R2.rstrip()
    cpcmd = "#cp ../%s_*.gz ./\nconda activate trinity2.9.1\n" % SRA
    conda = "conda activate trinity2.9.1\n"
    f_out.write(conda)

    f_out.write(cpcmd)
    # f_out.write(cmd1)
    # f_out.write(cmd2)
    dwnload.write(cmd1)
    dwnload.write(cmd2)
    # let trimmo these bad boys too in this script

    # q25 min len 57
    trimmo = "\n# QC trim the reads and remove vector\n#trimmomatic PE -summary trim_summary.txt  -threads 8 -phred33 %s_1.fastq.gz %s_2.fastq.gz %s_paired_1.fq crap1.fastq.gz %s_paired_2.fq crap2.fastq.gz ILLUMINACLIP:~/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:47\n" % (SRA, SRA, cell_line, cell_line)
    f_out.write(trimmo)
    delete_SRA_files = "#delete the raw reads\nrm *.fastq.gz\n"
    f_out.write(delete_SRA_files)
    # trinity de novo cmd
    tr_cmd = "\n# de nove RNAseq assembly: genome guided may not represent the cancer transcripts\n#Trinity --full_cleanup --seqType fq --left %s_paired_1.fq --right %s_paired_2.fq --CPU 8 --max_memory 250G --min_glue 3 --output %s_trinity\n" % (cell_line, cell_line, cell_line)
    f_out.write(tr_cmd)
    rename = "#cp %s_trinity/Trinity.fasta %s.trinity.fasta\n" % (cell_line, cell_line)
    f_out.write(rename)
    transrate = "\n#QC the assembly and trash stuff .. this is a bit gung ho!!!\n#transrate --assembly %s_trinity.Trinity.fasta --threads 8 --left %s_paired_1.fq --right %s_paired_2.fq\n" % (cell_line, cell_line, cell_line)
    f_out.write(transrate)
    f_out.write("\n# transdecoder to extract the cds\n")
    f_out.write("#TransDecoder.LongOrfs -t %s_trinity.Trinity.fasta\n" % cell_line)
    f_out.write("#diamond blastp --query ./*/longest_orfs.pep --db /gpfs1/scratch/bioinf/db/databases/uniprot.dmnd --threads 8 --max-target-seqs 1 --outfmt 6 -o blastp.outfmt6\n")
    f_out.write("#TransDecoder.Predict -t %s_trinity.Trinity.fasta --retain_blastp_hits blastp.outfmt6\n\n" % cell_line)

    
    f_out.write("# GENOME guided assembly - map the reads\n ")
    #STAR --runMode genomeGenerate --runThreadN 4 --limitGenomeGenerateRAM 259760745173 --genomeDir ./star_indicies --genomeFastaFiles Homo_sapiens.GRCh38.dna.toplevel.fa
    star_cmd = "STAR --genomeDir ../star_indicies/ --limitGenomeGenerateRAM 5554136874 --limitBAMsortRAM 5554136874 --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 7  --outFilterMultimapNmax 5 --outFileNamePrefix %s --readFilesIn  %s_paired_1.fq %s_paired_2.fq\n" %(cell_line, cell_line, cell_line)
    f_out.write(star_cmd)
    f_out.write("# GENOME guided assembly. Assemble from the mapped reads\n ")
    gg_cmd = "Trinity --genome_guided_bam %sAligned.sortedByCoord.out.bam --genome_guided_max_intron 10000 --max_memory 10G --CPU 12 --max_memory 100G --full_cleanup --output %s_GG_Trinity --genome_guided_min_coverage 5\n" % (cell_line, cell_line)
    f_out.write(gg_cmd)    
    



    # close this so we can open another file
    f_out.close()
dwnload.close()
    
