#!/usr/bin/env python
# title: Put SNP info in fasta description, keep original description
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
import errno
import logging
import logging.handlers
import sys
import difflib
from Bio import AlignIO # align those that fail
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices

# NOTE: converting NCBI human gene names to uniport:
# https://www.biostars.org/p/429062/

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
    sys.exit(1)

usage = """
VERSION = %s
to see the help: python detail_SNP_locations_in_fa_file.py -h

takes in the cds nt .fa for the human ref
and the reconstructed fa from a vcf. Plus the vcf used.

requires: Biopython, python 3.X. currently uses less than 2GB RAM
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
                          default="Homo_sapiens.GRCh38.cds.all.fa",
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

    optional.add_argument("--show", dest='show',
                          action="store",
                          default="NO",
                          type=str,
                          help="show the changes in those that fail")

    optional.add_argument("--align", dest='align',
                          action="store",
                          default="YES",
                          type=str,
                          help="show the alignment in those that fail")

    optional.add_argument("--remove_no_mut_genes", dest='remove_no_mut_genes',
                          action="store",
                          default="YES",
                          type=str,
                          help="do not output genes that are the same as reference")

    optional.add_argument("--min_len", dest='min_len',
                          action="store",
                          default=12,
                          type=int,
                          help="min len of AA seq to return. Some are as short as 2. We dont want these. ")

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
    if line.startswith("    # "):
        return False  # comment line
    if line.startswith("		p"):
        return False  # comment line
    return line.rstrip()


def split_line(line):
    """split line -  not used"""
    # was [1] for thapbi data
    # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
    return line.split()


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
                    # useful_info = INFO.split(";ANN=C|")[1]
                    useful_info = INFO.split(";ANN=")[1]
                    if (len(useful_info.split("|"))) < 4:
                        continue
                    try:
                        SNP_type = useful_info.split("|")[1]

                        SNP_impact = useful_info.split("|")[2]
                        gene_code = useful_info.split("|")[3]
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


def parse_NCBI_to_SP(infile, in_dict):
    """parse the file populate the dict"""
    with open(infile, 'r') as in_file:
        for line in in_file:
            if test_line(line):
                data = line.split("\t")
                ensem_gene_id = data[0]
                unipr_gn_id = ""
                unipr_gn_id = data[1]
                ensem_gene_id =  ensem_gene_id.split(",")
                for entry in ensem_gene_id:
                    unipr_gn_id = unipr_gn_id.rstrip()
                    in_dict[entry] = unipr_gn_id
    return in_dict


def parse_swiss_p_names(infile, in_dict):
    """parse the file populate the dict"""
    with open(infile, 'r') as in_file:
        for line in in_file:
            if test_line(line):
                data = line.split("|")
                sw_gn_id = data[1]
                in_dict[sw_gn_id] = line
    return in_dict


def parse_gene_info(gene_info, transcript_to_gene, in_dict):
    """column 9 of the human gff.
eg: ID=gene:ENSG00000187961;Name=KLHL17;biotype=protein_coding;
description=kelch like family member 17
[Source:HGNC Symbol%3BAcc:HGNC:24023];gene_id=ENSG00000187961;
logic_name=ensembl_havana_gene_homo_sapiens;version=14

    we want to extract
        ID=gene:ENSG00000187961
        Name=KLHL17

    returns a dictionary. gen name to actual name
    """

    data = gene_info.split(";")
    transcript = data[0].replace("ID=transcript:", "")
    gene = data[1].replace("Parent=gene:", "")
    name = data[2].replace("Name=", "")
    in_dict[transcript] = name
    in_dict[gene] = name
    transcript_to_gene[gene] = transcript
    transcript_to_gene[transcript] = gene
    return in_dict, transcript_to_gene



def parse_gff(gff):
    """function to parse GFF and produce a dictionary genes to names
    s..."""
    # iterate through gff
    gene_to_name = defaultdict(str)
    transcript_to_gene = defaultdict(str)
    f_in = open(gff)
    for line in f_in:
        if line.startswith("#"):
            continue
        if not line.strip():
                continue  #  if the last line is blank
        scaffold, prog, cds_type, start, stop, dot, direction, \
                  phase, gene_info = line.split("\t")
        if cds_type == "mRNA":
            gene_to_name, transcript_to_gene = parse_gene_info(gene_info,
                                                               transcript_to_gene,
                                                               gene_to_name)
    return gene_to_name, transcript_to_gene


def translate(sequences):
    """func to translate seq
    parse input file and save the AA in a dict
    returns
    [gene] = amino_acid_seq
    """
    gene_to_amino_acid = defaultdict(str)
    transcript_to_gene = defaultdict(str)
    gene_to_transcript = defaultdict(str)
    gene_names = []
    for seq_record in SeqIO.parse(sequences, "fasta"):
        gene_to_amino_acid[seq_record.id] = seq_record.seq.translate()
        gene_names.append(seq_record.id)
        elements = seq_record.description.split(" ")
        for i in elements:
            if "gene:" in i:
                gene = i.split("gene:")[1]
                transcript_to_gene[seq_record.id] = gene
                transcript_to_gene[seq_record.id] = gene.split(".")[0]
                gene_to_transcript[gene] = seq_record.id
                gene_to_transcript[gene] = seq_record.id.split(".")[0]

    return gene_to_amino_acid, gene_names, transcript_to_gene, gene_to_transcript


def return_codon(seq, pos):
    """func to return the codon of the nt snp of interest
    takes in :
    sequnce of interest. The post (int) of the nt SNP.
    return the codon, positionally aware of where the SNP was.
    """
    if pos %3 == 0:
        # this SNP is the first position in the codon
        return seq[pos:(pos+3)]
    if pos %3 == 1:
        # this SNP is the 2nd position in the codon
        return seq[(pos -1):(pos+2)]
    if pos %3 == 2:
        # this SNP is the 3nd/last position in the codon
        return seq[(pos -2):(pos+1)]



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

    # check files exist or not
    for user_file in file_list:
        if not os.path.isfile(user_file):
           logger.warning("file not found: %s", user_file)

    # create a dictionary of NCBI gene to uniprot names
    NCBI_to_prot_file = os.path.join("db", "transcript_names.txt")
    NCBI_to_prot_dict = defaultdict(str)

    NCBI_to_prot_dict = parse_NCBI_to_SP(NCBI_to_prot_file,
                                         NCBI_to_prot_dict)
    prot_file = os.path.join("db", "swiss_prot_ids")
    description_to_swiss_id = defaultdict(str)


    description_to_swiss_id = parse_swiss_p_names(prot_file,
                                                  description_to_swiss_id)

    logger.info("swiss prot names loaded")
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
    logger.info((dict(list(gene_to_snp_type.items())[:1])))


    logger.info("parsing gff file")
    gene_to_name, transcript_to_gene_name = parse_gff(args.gff)

    # collect the REF protein seq
    REF_gene_to_amino_acid, wt_gene_names, \
            transcript_to_gene, gene_to_transcript = translate(args.wt)
    # collect the ALT protein seq
    ALT_gene_to_amino_acid, alt_gene_names, \
            transcript_to_gene, gene_to_transcript = translate(args.cancer_fa)

    # iterate through the genes and compare what bcftools did
    # bcftools will not output all variants, so using vcf here
    # may not help!!
    # open the outfile
    f_out = open(args.out, "w")
    logger.info("comparing REF amino acids to ALT amino acids")
    # dict to capture to relavant info
    gene_description_location = defaultdict(set)
    seen_gene = set()
    seen_gene2 = set()
    logger.info("sets set up")

    # lets write the alignmne to a file, generated in next loop
    align_out = args.cancer_fa.split(".fa")[0] + "_vs_WT.align.fasta"
    f_align = open(align_out, "w")
    logger.info("outfiles set up")

    # debug file to extract the nucleotide variant responsible.
    debug = open("debug.txt", "w")

    # iterate through
    logger.info("now to iterate through")
    written_set = set([])
    failed_set =  set([])
    for wt_record in SeqIO.parse(args.wt, "fasta"):
        
        gene = wt_record.id

        # get the nucleotide seqs
        if not gene in wt_nt:
            continue
        wt_seq = wt_nt[gene]
        wt_sequence = str(wt_seq.seq)

        # for the alt
        alt_seq = alt_nt[gene]
        alt_sequence = str(alt_seq.seq)

        # get the AA seq
        wt_seq = REF_gene_to_amino_acid[gene]
        # get the mutant seq
        alt_seq = ALT_gene_to_amino_acid[gene]

        if wt_seq.upper() != alt_seq.upper():
            # use alignment to find differences
            blosum62 = substitution_matrices.load("BLOSUM62")
            # see: https://github.com/biopython/biopython/issues/3634
            # ds not xx
            # blosum62, -10, -0.5
            alignments = pairwise2.align.globalds(wt_seq, alt_seq, blosum62, -10, -1)

            # logger.info(format_alignment(*alignments[0]))
            differences = []
            pos = 0
            wt_gap = 0
            alt_gap = 0

            # lets write the alignent to a file
            out_align = ">WT (top) %s _vs_ cancer %s\n" %(gene, gene)
            f_align.write(out_align)
            f_align.write(format_alignment(*alignments[0], full_sequences=True))

            for wt, alt in zip(alignments[0][0],alignments[0][1]):
                pos += 1
                if wt == alt:
                    continue
                else:
                    if wt == "-":
                        wt_gap += 1
                    if alt == "-":
                        alt_gap += 1
                    data = "%d\t%s\t%s\t%d\t%d" % (pos, wt, alt,
                                                   wt_gap, alt_gap)
                    #if wt != "-" or alt != "-":
                        #muta_info = ("WT\t%s\tALT\t%s" % (wt, alt))
                        #print(muta_info)
                    differences.append(data)
                #print(differences)
            protien_sub = ""
            for change in differences:
                pos, wt, alt, wt_gap, alt_gap = change.split("\t")
                count = int(pos) - 1
                out_data =  ""
                try:
                    count = count -1 # for indexing
                    wt_seq[count] # for some reason sometimes this is beyond len(gene)

                    #if wt_seq[count] == alt_seq[count]: # sometimes it is lying to me!!
                        #continue
                    if wt == "-":
                        wt = "del"
                    if alt == "-":
                        alt = "del"

                    protien_s = "p_%s%d%s " % (wt, count,
                                              alt)
                    protien_sub =  protien_sub + protien_s
                    nt_count = (count +1)*3
                    nucleo_sub = "c_%s%d%s" % (wt_sequence[nt_count],
                                               nt_count,
                                               alt_sequence[nt_count])
                    # nucleo_sub = "" # this doesnt work due to "-" in alignments

                    # seq of interest would be the alt, right?
                    seq_of_int = alt_seq[count - 4: count + 4]
                    original_seq = wt_seq[count - 4: count + 4]

                    #################################################################
                    # debug to find correct nt seq
                    wt_codon = return_codon(wt_sequence, nt_count)
                    alt_codon = return_codon(alt_sequence, nt_count)

                    part1 = "gene: %s\tWT: %s\tALT: %s\t nt_pos:\t%d\t" % (gene,
                                                            wt, alt, nt_count)
                    part2 = "nt:\t%s\torignal: %s\t alt:\t%s\t" % (wt_sequence[nt_count],
                                                            original_seq, seq_of_int)
                    part3 = "original_codon:\t%s alt_codon:\t%s\n" % (wt_codon, alt_codon)

                    debug_out = part1 + part2 + "\n" # + part3

                    debug.write(debug_out)

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
                    this_fails = "warning %s fails" % gene #  usually here - if not wt_seq[count]:
                    logger.info(this_fails)
                    failed_set.add(wt_record.id)
                alt_seq_record_AA = ALT_gene_to_amino_acid[gene]
                # get the swiss prot gene name from the NCBI transcript name
                swiss_prt = NCBI_to_prot_dict[wt_record.id.split(".")[0]]
                # get the swiss prot info:
                swiss_prto_info = description_to_swiss_id[swiss_prt]
                if swiss_prto_info != "":
                    wt_record.id = swiss_prto_info + " | "
                #remove the extra info as requested.
                temp = str(wt_record.description)
                try:
                    temp = temp.split(" description:")[1]
                    temp = temp.split("[Source:")[0]
                except:
                    temp = ""
                wt_record.description = " | " + temp  + " | " + out_data
                # swap the AA seq for the variant type. 
                wt_record.seq =  alt_seq_record_AA
            if len(wt_record.seq) > args.min_len:
                if gene not in written_set:
                    print("writing .. %s" % gene)
                    wt_record.id = wt_record.id.replace("|   |  |", "|")
                    temp = wt_record.description
                    wt_record.description = temp.replace("|   |  |", "|")
                    SeqIO.write(wt_record, f_out, "fasta")
                    written_set.add(gene)


    logger.info("writting out the failed alignment seqs")
    for wt_record.id in failed_set:
        gene = wt_record.id
        wt_record = wt_nt[gene]
        alt_seq_record_AA = ALT_gene_to_amino_acid[gene]
        # get the swiss prot gene name from the NCBI transcript name
        swiss_prt = NCBI_to_prot_dict[wt_record.id.split(".")[0]]
        # get the swiss prot info:
        swiss_prto_info = description_to_swiss_id[swiss_prt]
        if swiss_prto_info != "":
            wt_record.id = swiss_prto_info + " | ALIGNMENT_FAIL_COMPLEX_CHANGES |"
        #remove the extra info as requested.
        temp = str(wt_record.description)
        try:
            temp = temp.split(" description:")[1]
            temp = temp.split("[Source:")[0]
        except:
            temp = ""
        # swap the AA seq for the variant type. 
        wt_record.seq =  alt_seq_record_AA
        if len(wt_record.seq) > args.min_len:
            if gene not in written_set:
                print("writing .. %s" % gene)
                wt_record.id = wt_record.id.replace("|   |  |", "|")
                temp = wt_record.description
                wt_record.description = temp.replace("|   |  |", "|")
                SeqIO.write(wt_record, f_out, "fasta")
                written_set.add(gene)
    # close fasta out, close align out
    f_out.close()
    f_align.close()
    debug.close()
    logger.info("finished")






