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
    read1 = R1.split("/")[-1]
    read2 = R2.split("/")[-1]
    SRA = R1.split("/")[-1]
    SRA = SRA.split("_1")[0]
    cell_line = cell_line.rstrip()
    outsh = "%s_seq2hla_raw.sh" % cell_line.rstrip()
    #print("fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files %s\n" %SRA)

    f_out = open(outsh, "w")
    # 6 cores, so only 1 runs on a "node*" at a time. 
    print("\nsbatch %s_seq2hla_raw.sh\n" % cell_line)
    f_out.write("#!/bin/bash -l\n")
    f_out.write("#SBATCH -J %s_RNAseq\n" % cell_line)
    f_out.write("#SBATCH --tasks-per-node=16\n")
    f_out.write("#SBATCH -p bigmem\n")
    f_out.write("#SBATCH --mem=60GB\n")
    f_out.write("cd /gpfs1/scratch/bioinf/pjt6/cancer_genomes/seq2hla\n")

    mkdir = "# make a folder of the cell name\nmkdir %s_Raw\n" % cell_line
    f_out.write(mkdir)
    cd_into = "cd ./%s_Raw\n" % cell_line
    f_out.write(cd_into)

    conda = "conda activate seq2HLA\n"
    f_out.write(conda)

    # let trimmo these bad boys too in this script


    f_out.write("python /home/pjt6/scratch/cancer_genomes/seq2hla/seq2HLA/seq2HLA.py -1 /home/pjt6/scratch/cancer_genomes/fq/%s -2 /home/pjt6/scratch/cancer_genomes/fq/%s -r %s_raw -p 16 \n" % (read1.split(".gz")[0],read2.split(".gz")[0], cell_line))




    # close this so we can open another file
    f_out.close()
    
