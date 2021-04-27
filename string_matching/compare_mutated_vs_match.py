# title: Script to macth coloumn for colomn mutated hits.
# this does not implemen fuzzymatching.
# author. P. Thorpe April 2021

import sys
import re
import time
import os
import pandas as pd
import numpy as np
from fuzzysearch import find_near_matches
from collections import defaultdict
from collections import OrderedDict
from optparse import OptionParser
import errno
import logging
import logging.handlers


def split_up_mut_pep(peptide, kmer=4):
    """function to split up the mutant protein.
    at the moment into 5 chunks based on kmer length"""
    for pep in range(0,len(peptide)):
        beg = peptide[0:kmer]
        end1 = peptide[kmer-1:(kmer+kmer-1)]
        end2 = peptide[(len(peptide)-kmer):(len(peptide))]
        mid_1 = peptide[1:5]
        mid_2 = peptide[2:6]
        combinations = [beg, end1, end2, mid_1, mid_2]
    return combinations


def parse_excel(infile, out_file):
    """use pandas to open excel.
    Need coloumns called:
     Mutated Sequence  and Peptide Sequence
    """
    in_data = pd.read_excel(infile)
    out_df = pd.DataFrame()
    hitdf = pd.DataFrame()
    #print(in_data)
    sequences = in_data['Peptide']
    mutated = in_data['Mutated']
    line = 1
    hits = OrderedDict()
    for i in range(0, len(in_data.index)):
        
        mutated_seq = mutated[line -1]
        mutated_seq = mutated_seq.strip()
        seq =  sequences[line -1]
        seq = seq.strip()
        combinations = split_up_mut_pep(mutated_seq, kmer=4)
        last_frag = ""
        for pep_frag in combinations:
            if pep_frag == last_frag:
                continue
            last_frag = pep_frag
            if pep_frag in seq:
                print(line -1, "\t", mutated_seq, "\t", seq)
                MM_s = find_near_matches(pep_frag, seq,
                                         max_deletions=0,
                                         max_insertions=0,
                                         max_substitutions=0,
                                         max_l_dist=1)
                out_df = out_df.append(in_data.iloc[[i]])
                if MM_s != []:
                    name = str(line-1) + "_" + pep_frag + "_" + seq
                    hits[name] = MM_s      
        line += 1
    out_df.to_excel(out_file)
    return hits


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.2")
    sys.exit(0)

if sys.version_info[:1] != (3,):
    # e.g. sys.version_info(major=3, minor=6, micro=7,
    # releaselevel='final', serial=0)
    # break the program
    print ("currently using:", sys.version_info,
           "  version of python")
    raise ImportError("Python 3.x is now required for this.py")
    print ("did you activate the virtual environment?")
    sys.exit(1)



usage = """Use as follows:
$ python compare_mutated_vs_match.py -i input_excel.xlxs -o output_excel.xlxs

input file can have extra coloumns but must have:

coloumns called: Mutated Sequence  and Peptide Sequence
"""

parser = OptionParser(usage=usage)

parser.add_option("-i",
                  dest="in_excel",
                  default="fuzzy_test.xlsx",
                  help="this is the input excel. Need coloumns called." +
                  " Mutated Sequence  and Peptide Sequence.")

parser.add_option("-o", "--out",
                  dest="outfile",
                  default="test.xlsx",
                  help="output data in excel  " +
                  " what name do you want?")


# get the user options. TODO. arg parser instead
(options, args) = parser.parse_args()
in_excel = options.in_excel
outfile = options.outfile

logfile = in_excel.split(".xa")[0] + ".log"

# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('compare_mutated_vs_match.py: %s'
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

    if not os.path.isfile(in_excel):
        print("file not found: %s" % user_file)
        os._exit(0)
    outfmt = "looking at: %s" % in_excel
    logger.info(outfmt)
    hits = parse_excel(in_excel, outfile)
    outfile = outfile.split(".x")[0] + "hit_locations.txt"
    f_out = open(outfile, "w")
    f_out.write("#line\tMUT_PEP\tSequence\tMatched\n")
    for key, items in hits.items():
        print(key, items)
        line = key.split("_")[0]
        mut = key.split("_")[1]
        sequences =  key.split("_")[2]
        outdata = "%s\t%s\t%s\t%s\n" % (line, mut, sequences, str(items))
        f_out.write(outdata)
    f_out.close()
