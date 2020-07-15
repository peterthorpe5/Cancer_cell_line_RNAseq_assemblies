#!/usr/bin/env python
# title: compare MS results for WT and mutant
# (c) The James Hutton Institute 2020
# Author: Peter Thorpe
import sys
import os
import argparse
from collections import defaultdict


VERSION = "summerise results: v0.01"
if "--version" in sys.argv:
    print(VERSION)
    sys.exit(1)

usage = """
%s

# to see the help


 python .py -h
 
""" % VERSION

if "--help" or "-h" in sys.argv:
    print(usage)


def get_args():
    parser = argparse.ArgumentParser(description="Pipeline: summerise results " +
                                     "in a folder ",
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("compare")[0]
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--wt", dest='wt_ms',
                          action="store",
                          default="HCC 1395 Trial_PEAK DB on de novo 2_UP56400_protein-peptides.csv",
                          type=str,
                          help="WT MS data e.g. HCC 1395 Trial_PEAK DB on de novo 2_UP56400_protein-peptides.csv ")


    optional.add_argument("--snp", dest='mt_ms',
                          action="store",
                          default="HCC 1395 Trial_Re Run PEAK DB on De novo only_SNP customised database_protein-peptides.csv",
                          type=str,
                          help="snp ms data: HCC 1395 Trial_Re Run PEAK DB on De novo only_SNP customised database_protein-peptides")

    optional.add_argument("--wt_fa", dest='wt_fa',
                          action="store",
                          default="Homo_sapiens.GRCh38.cds.all.fa",
                          type=str,
                          help="WT fa: Homo_sapiens.GRCh38.cds.all.fa")

    optional.add_argument("--mt_fa", dest='mt_fa',
                          action="store",
                          default="HCC_1395_alt_from_SNPs.fasta",
                          type=str,
                          help="MT fa: HCC_1395_alt_from_SNPs.fasta")

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
    return line



def parse_file(infile, gene_to_prot, prot_start_stop):
    """funct to parse the MS data file"""
    with open(infile, "r") as fh:
        for line in fh:
            line = split_line(line)
            if not test_line(line):
                continue
            if line.startswith("Protein"): # 1 st line
                pass
            Protein_Group, Protein_ID, Protein_Accession, Peptide, Unique,\
                           ten_10lgP, Mass, Length, ppm, m_z, z, RT, \
                           Area_Sample_1, Fraction, Scan, Source_File, \
                           Feature, Feature_Sample_1, \
                           Start, End, PTM, AScore, \
                           Found_By = line.split(",")
            gene_to_prot[Protein_Accession].add(Peptide)
            gene_plus_prot = "%s_%s" % (Protein_Accession, Peptide)
            start_stop = "%s_%s" % (Start, End)
            prot_start_stop[gene_plus_prot] = start_stop
    return gene_to_prot, prot_start_stop            



args, FILE_DIRECTORY = get_args()
# Run as script
if __name__ == '__main__':
    # collect data for the WT 
    wt_gene_to_prot = defaultdict(set)
    wt_prot_start_stop = defaultdict(str)
    # gene_plus_prot = "%s_%s" % (Protein_Accession, Peptide)
    wt_gene_to_prot, wt_prot_start_stop = parse_file(args.wt_ms,
                                                     wt_gene_to_prot,
                                                     wt_prot_start_stop)

        # collect data for the mt 
    mt_gene_to_prot = defaultdict(set)
    mt_prot_start_stop = defaultdict(str)
    # gene_plus_prot = "%s_%s" % (Protein_Accession, Peptide)
    mt_gene_to_prot, mt_prot_start_stop= parse_file(args.mt_ms,
                                                    mt_gene_to_prot,
                                                    mt_prot_start_stop)
    for gene, peptides in mt_gene_to_prot.items():
        print(gene, "peptides", peptides)
    
 
     



