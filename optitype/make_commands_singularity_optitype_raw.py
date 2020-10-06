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
    print("\nsingularity run /storage/home/users/pjt6/docker/optitype/optitype_latest.sif --rna --input %s_1.fastq.gz %s_2.fastq.gz -o %s_raw_optitype -v\n" % (SRA,
                                                                                                                                                                    SRA,
                                                                                                                                                                    cell_line))
