#!/usr/bin/env python
# title: Put SNP info in fasta description
# (c) The James Hutton Institute 2021
# Author: Peter Thorpe

# import
import sys
import os
import argparse
from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO
import time
import os
import errno
import logging
import logging.handlers
import sys
import difflib


VERSION = "find SNP locations: v0.01"
if "-v" in sys.argv or "--version" in sys.argv:
    print(VERSION)
    sys.exit(0)

if sys.version_info[:1] != (3,):
    # e.g. sys.version_info(major=3, minor=6, micro=7,
    # releaselevel='final', serial=0)
    # break the program
    print ("currently using:", sys.version_info,
           "  version of python")
    raise ImportError("Python 3.x is now required for this .py")
    print ("did you activate the virtual environment?")
    print ("this is to deal with module imports")
    sys.exit(1)

usage = """
VERSION = %s
to see the help: python detail_SNP_locations_in_fa_file.py -h

takes in the cds nt .fa for the human ref
and the reconstructed fa from a vcf. Plus the vcf used.

requires: Biopython, python 3.X. currently uses less than 1GB RAM
""" % VERSION

if "--help" or "-h" in sys.argv:
    print(usage)


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: summerise results " +
                                     "in a folder ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("compare")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--wt", dest='wt',
                          action="store",
                          default="Homo_sapiens.GRCh38.cds.all.fa2",
                          type=str,
                          help="the ref human proteins ")

    optional.add_argument("--vcf", dest='vcf',
                          action="store",
                          default="MDAMB231_filtered.q30.mindp3.vcf.recode2",
                          type=str,
                          help="vcf file with the snps detailed. Has to be against the genes only" +
                          " this will not do any filtering")
    
    optional.add_argument("--snp_eff_vcf", dest='snp_vcf',
                          action="store",
                          default="MDAMB231_Snpeff.vcf",
                          type=str,
                          help="vcf file with the snps detailed. Has to be against the genes only" +
                          " this will not do any filtering")

    optional.add_argument("--gff", dest='gff',
                          action="store",
                          default="Homo_sapiens.GRCh38.99.gff3",
                          type=str,
                          help="Homo_sapiens.GRCh38.99.gff3 used to convert gene name " +
                          " why this is so difficult I will never know")

    optional.add_argument("--cancer_fa", dest='cancer_fa',
                          action="store",
                          default="MDAMB231_H1_SNPs_reconstruct.q30.mindp3.fasta2",
                          type=str,
                          help="the reconstructed fa using the ref and vcf using bcftools")

    optional.add_argument("-o", "--out", dest='out',
                          action="store",
                          default="results-summary.txt",
                          type=str,
                          help="the tab file to fill in. Input file")

    optional.add_argument("-h", "--help",
                          action="help",
                          default=argparse.SUPPRESS,
                          help="Displays this help message"
                          " type --version for version")

    optional.add_argument('--version',
                          action='version',
                          version="%s: populate.py " + VERSION)
    args = parser.parse_args()
    return args, file_directory


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    if line.startswith("    # "):  # swarm result file
        return False  # comment line
    if line.startswith("		p"):
        return False  # comment line
    return line.rstrip()


def split_line(line):
    """split line"""
    # was [1] for thapbi data
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO

    return line.split


def parse_snpeff_vcf(vcf):
    """func to parse vcf file
    in: vcf file
    out: dictionary
        [gene_POS] = Gene \t POS \t REF \t ALT
    """
    gene_to_snp_type = defaultdict(list)
    with open(vcf, 'r') as v_file:
        for line in v_file:
            if test_line(line):
                chromo, POS, ID, REF, ALT, \
                       QUAL, FILTER, INFO, \
                       FORMAT, unknown = line.split()
                # print(INFO.split(";ANN=C|"))
                try:
                    useful_info = INFO.split(";ANN=C|")[1]
                    if (len(useful_info.split("|"))) < 4:
                        continue
                    try:
                        SNP_type = useful_info.split("|")[0]
                        
                        SNP_impact = useful_info.split("|")[1]
                        gene_code = useful_info.split("|")[2]
                        gene_pos = gene_code + "_" + POS
                        data = "%s\t%s\t%s\t%s" % (gene_code,
                                                   POS, SNP_type,
                                                   SNP_impact)
                        gene_to_snp_type[gene_code].append(data)
                        # BREAKING HERE!! GENE?
                        # print(gene_code, SNP_type)
                    except:
                        continue
                except:
                    continue
    return gene_to_snp_type


def parse_vcf(vcf):
    """func to parse vcf file
    in: vcf file
    out: dictionary
        [gene_POS] = Gene \t POS \t REF \t ALT
    """
    gene_snp_location = defaultdict(str)
    with open(vcf, 'r') as v_file:
        for line in v_file:
            if test_line(line):
                GENE, POS, ID, REF, ALT, \
                       QUAL, FILTER, INFO, \
                       FORMAT, unknown = line.split()
                data = "%s\t%s\t%s\t%s" % (GENE, POS, REF, ALT)
                gene_pos = GENE + "_" + POS
                gene_snp_location[gene_pos] = data
    return gene_snp_location


def parse_gene_info(gene_info, in_dict):
    """column 9 of the human gff.
eg: ID=gene:ENSG00000187961;Name=KLHL17;biotype=protein_coding;description=kelch like family member 17 [Source:HGNC Symbol%3BAcc:HGNC:24023];gene_id=ENSG00000187961;logic_name=ensembl_havana_gene_homo_sapiens;version=14

    we want to extract ID=gene:ENSG00000187961
    Name=KLHL17
    """
    data = gene_info.split(";")
    transcript = data[0].replace("ID=transcript:", "")
    gene = data[1].replace("Parent=gene:", "")
    name = data[2].replace("Name=", "")
    in_dict[transcript] = name
    in_dict[gene] = name
    return in_dict


def parse_gff(gff):
    """function to parse GFF and produce a dictionary genes to names
    s..."""
    # iterate through gff
    gene_to_name = defaultdict(str)
    f_in = open(gff)
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        scaffold, prog, cds_type, start, stop, dot, direction, \
                  phase, gene_info = line.split("\t")
        if cds_type == "mRNA":
            gene_to_name = parse_gene_info(gene_info, gene_to_name)
    return gene_to_name


def translate(sequences):
    """func to translate seq
    parse input file and save the AA in a dict
    returns
    [gene] = amino_acid_seq
    """
    gene_to_amino_acid = defaultdict(str)
    gene_names = []
    for seq_record in SeqIO.parse(sequences, "fasta"):
        gene_to_amino_acid[seq_record.id] = seq_record.seq.translate()
        gene_names.append(seq_record.id)
    return gene_to_amino_acid, gene_names


# get the agrs
args, FILE_DIRECTORY = get_args()
# set up a logging file
logfile = args.cancer_fa.split(".fa")[0] + ".WARNINGS.log"


# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('detail_SNPs.py: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    # Report input arguments
    logger.info(sys.version_info)
    logger.info("Command-line: %s", ' '.join(sys.argv))
    logger.info("Starting: %s", time.asctime())
    
    # check the files exist
    file_list = [args.wt, args.cancer_fa,
                 args.vcf, args. snp_vcf,
                 args.gff]
    
    for user_file in file_list:
        if not os.path.isfile(user_file):
           logger.warning("file not found: %s", user_file)
           os._exit(0)
    # collect data for the WT .fa
   
    logger.info("Indexing... %s", args.wt)
    wt_nt =  SeqIO.index(args.wt, "fasta")
    # collect data for the reconstructed .fa
    logger.info("Indexing... %s", args.cancer_fa)
    alt_nt =  SeqIO.index(args.cancer_fa, "fasta")

    # [gene_POS] = Gene \t POS \t REF \t ALT
    logger.info("Indexing... %s", args.vcf)
    gene_snp_location = parse_vcf(args.vcf)
    logger.info("Indexing... %s", args.snp_vcf)
    gene_to_snp_type = parse_snpeff_vcf(args.snp_vcf)
   
    # print(gene_to_snp_type)
    # quick look at the dict
    logger.info("SNP type quick look")
    logger.info((dict(list(gene_to_snp_type.items())[:3])))

    
    logger.info("parsing gff file")
    gene_to_name = parse_gff(args.gff)

    # collect the REF protein seq
    REF_gene_to_amino_acid, wt_gene_names = translate(args.wt)
    # collect the ALT protein seq
    ALT_gene_to_amino_acid, alt_gene_names = translate(args.cancer_fa)

    # iterate through the genes and compare what bcftools did
    # bcftools will not output all variants, so using vcf here
    # may not help!!
    # open the outfile
    f_out = open(args.out, "w")
    logger.info("comparing REF amino acids to ALT amino acids")
    # dict to capture to relavant info
    gene_description_location = defaultdict(set)
    seen_gene = set()

    # iterate through 
    for gene_pos, info in gene_snp_location.items():
        gene, POS, REF, ALT = info.split()
        if gene in seen_gene:
            continue
        seen_gene.add(gene)
        
        # get the nucleotide seqs
        wt_seq = wt_nt[gene]
        wt_sequence = str(wt_seq.seq)
        
        # for the alt
        alt_seq = alt_nt[gene]
        alt_sequence = str(alt_seq.seq)

        # get the AA seq
        wt_seq = REF_gene_to_amino_acid[gene]
        # get the mutant seq
        alt_seq = ALT_gene_to_amino_acid[gene]
        if wt_seq.upper() == alt_seq.upper():
            #alt_seq_record = ALT_gene_to_amino_acid[gene]
            #SeqIO.write(alt_seq_record, f_out, "fasta")
            continue
        for count, s in enumerate(difflib.ndiff(wt_seq, alt_seq)): #  this takes too long
            if s[0] == ' ':
                continue
            if not s[0]:
                continue
            elif s[0] == '-':
                #print("%s\talternative\t%s\tfrom position\t%d" % (gene, s[-1], count))
                try:
                    wt_seq[count] # for some reason sometimes this is beyond len(gene)

                    if wt_seq[count] == alt_seq[count]: # sometimes it is lying to me!!
                        continue
                    protien_sub = "p_%s%d%s" % (wt_seq[count], count,
                                                alt_seq[count])
                    nt_count = (count +1)*3
                    nucleo_sub = "c_%s%d%s" % (wt_sequence[nt_count],
                                               nt_count,
                                               alt_sequence[nt_count])
                    seq_of_int = wt_seq[count - 4: count + 4]
                    # print(protien_sub)
                    # need to get the stupid other names
                    # print(gene.split(".")[0])
                    name = gene_to_name[gene.split(".")[0]]
                    # snp eff name do not have - in them
                    name =  name.split("-")[0]
                    # print(name)
                    # get the SNPeff stuff
                    snp_eff_data = gene_to_snp_type[name]
                    # print("gene", gene, "name", name)
                    #print(snp_eff_data)
                    temp = ""
                    for i in snp_eff_data:
                        temp = temp + " " + i
                        # print(i)
                    out_data = "\t%s|%s|%s|%s|" % (nucleo_sub, protien_sub,
                                                   seq_of_int, temp)
                    gene_description_location[gene].add(out_data)
                except:
                    this_fails = "warning" #  usually here - if not wt_seq[count]:
                    continue

    logger.info("finished comparing. Will now output the data")
    for gene in alt_gene_names:
        alt_seq_record = alt_nt[gene]
        alt_seq_record_AA = ALT_gene_to_amino_acid[gene]
        alt_seq_record.seq = alt_seq_record_AA
        # print(alt_seq_record)
        SNP_info = gene_description_location[gene]
        temp = ""
        for entry in SNP_info:
            temp =  temp + entry
        alt_seq_record.description = alt_seq_record.description  + temp
        SeqIO.write(alt_seq_record, f_out, "fasta")
    f_out.close()
    logger.info("finished")






