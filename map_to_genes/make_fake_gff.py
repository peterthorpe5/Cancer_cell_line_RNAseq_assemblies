#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO



def seq_getter(filename1, outfile):
    "opens up a fasta files and prints all names of sequences to a file."

    f_out = open(outfile, "w")
    f_out.write("##gff-version 3\n")
    f = open("aa.fa", "w")
    for seq_record in SeqIO.parse(filename1, "fasta"):
        outdata = "%s\tGenBank\tgene\t1\t%d\t.\t+\t.\tID=%s\n" % (seq_record.id,
                                                                  len(seq_record.seq),
                                                                  seq_record.description)
        
        f_out.write(outdata)
        seq_record.seq = seq_record.seq.translate()
        seq_record.id = seq_record.description
        SeqIO.write(seq_record, f, "fasta")
    f_out.close()
    f.close()


usage = """Use as follows:
$ python re_format_gff.py --gff augustus.gff -s src=E (default) -o agustus_reformatted
script to reformt gff coloumn 9 as the ; formatting brakes some tools. 
"""

parser = OptionParser(usage=usage)


parser.add_option("-i", "--in", dest="filename1", default=None,
                  help="in fasta'",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="outfile", default=None,
                  help="outfile to the output filenames")


(options, args) = parser.parse_args()


filename1 = options.filename1
outfile = options.outfile

seq_getter(filename1, outfile)
