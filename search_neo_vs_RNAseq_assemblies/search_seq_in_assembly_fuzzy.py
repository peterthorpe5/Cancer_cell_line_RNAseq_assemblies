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
from Bio.Seq import Seq
from Bio import SeqIO
import time
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_write_out6frame(infasta, outfasta):
    """func to parse a fasta nt file and write
    out the 6 frame translation of this"""

    f_out = open(outfasta, 'w')
    for rec in SeqIO.parse(infasta, 'fasta'):
        for frame in range(0,3):
            fram_des = "frame %d" % frame
            SeqIO.write(SeqRecord(seq=rec[frame:].seq.translate(),
                                  id=rec.id,
                                  description=fram_des),
                        f_out, 'fasta')
        rec.seq = rec.seq[::-1]
        for frame in range(0,3):
            fram_des = "frame - %d" % frame
            SeqIO.write(SeqRecord(seq=rec[frame:].seq.translate(),
                                  id=rec.id,
                                  description=fram_des),
                        f_out, 'fasta')
    f_out.close()

def split_seq_into_peptide_size_chunks(string, slice_len=10):
    size = len(string)
    chunks_list = []
    start = 0
    for i in range(size):
        end = start + slice_len -1
        slice = string[start:end]
        print("start", start, "end", end)
        print(slice)
        start =  start + 1
        chunks_list.append(slice)
   
def fa_out6frame(fasta, cell_neoep, logger):
    """func to parse the 6 frame trans
     the name  """
    collect_hits = set([])
    cell_neo = dict()
    neo_set = set([])
    data = []
    filename = fasta
    out_set = set([])
    for entry in cell_neoep:
        line, neo = entry.split("\t")
        cell_neo[neo] = line
        neo_set.add(neo)
    for seq_record in SeqIO.parse(fasta, "fasta"):
        sequence = seq_record.seq
        for peptide in neo_set:
            peptide_len = len(peptide)
            # sequence_split_list = split_seq_into_peptide_size_chunks(seq_record.seq, peptide_len)
            cell_line = cell_neo[peptide]
            temp = peptide.rstrip()
            test_pep = Seq(peptide)
            if peptide in seq_record.seq:
                result = "\t".join(["line: ", cell_line, "peptide",
                            peptide, "hit in", filename, 
                            seq_record.id, seq_record.description, 
                            str(data)])
                if result not in out_set:
                    logger.info(result)
                    out_set.add(result)
            #for fragment in sequence_split_list:
                #if len(fragment) < peptide_len:
                    #continue
#            data = find_near_matches(seq_record.seq, test_pep,
#            max_deletions=1, max_insertions=1, max_substitutions=1,
#            max_l_dist=0)
#            
#            if data != []:
#                print("hehre cell line = ", cell_line)
#                result = "\t".join(["line: ", cell_line, "peptide",
#                            peptide, "hit in", str(fragment), filename, 
#                            seq_record.id, seq_record.description, 
#                            str(data)])
#                logger.info(result)

                #collect_hits.add(result)
         

cell_neoep = """HCC 1806	RQVTSSGVSY
HCC 1806	RQYLQEVGY
HCC 1806	ESGPSIVHR
HCC 1954	YPDRIMNTF
HCC 1954	KYILSNANLF
HCC 1954	KAYHEQLTV
HCC 1954	FEQEMATAA
HCC 1954	EYPDRIMNTF
MDA-MB-231	KYLDEDTIYHL
MDA-MB-415	DEFNVQVL
MDA-MB-415	DNFLMGIGR
MDA-MB-415	ENIIFEEY
MDA-MB-453	LLDSSQKNLY""".split("\n")


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
$ python search_seq_in_assembly.py -i fa -o outfile

make 6 frame translation of the assemlby
then serahces for the peptides, allowing fuzzyness. for now
"""

parser = OptionParser(usage=usage)

parser.add_option("-i",
                  dest="in_fasta",
                  default=None,
                  help="de novo assembly" +
                  " Mutated Sequence  and Peptide Sequence.")

parser.add_option("-o", "--out",
                  dest="outfile",
                  default=None,
                  help="6 frame translation")


# get the user options. TODO. arg parser instead
(options, args) = parser.parse_args()
in_fasta = options.in_fasta
outfile = options.outfile

logfile = in_fasta.split(".xa")[0] + ".log"

# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('search_seq_in_assembly.py: %s'
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

    if not os.path.isfile(in_fasta):
        print("file not found: %s" % user_file)
        os._exit(0)
    outfmt = "looking at: %s" % in_fasta
    logger.info(outfmt)
    logger.info("translating sequences")
    # parse_write_out6frame(in_fasta, outfile)
    logger.info("line to peptide info")
    # logger.info(cell_neoep)
    logger.info("searching sequences")
    fa_out6frame(outfile, cell_neoep, logger)

cmd = """
python search_seq_in_assembly_fuzzy.py -i HCC_1806_trinity.Trinity.fasta -o HCC_1806_trinity_6frame
python search_seq_in_assembly_fuzzy.py -i MDAMB231_trinity.Trinity.fasta -o MDAMB231_trinity_6frame
python search_seq_in_assembly_fuzzy.py -i MDAMB453_trinity.Trinity.fasta -o MDAMB453_trinity_6frame
python search_seq_in_assembly_fuzzy.py -i HCC1954_trinity.Trinity.fasta -o HCC1954_trinity_6frame
python search_seq_in_assembly_fuzzy.py -i MDAMB415_trinity.Trinity.fasta -o MDAMB415_trinity_6frame


python search_seq_in_assembly_fuzzy.py -i HCC_1806_trinity.Trinity.fasta.transdecoder.cds -o HCC_1806_trinity.Trinity.fasta.transdecoder.pep 

python search_seq_in_assembly_fuzzy.py -i MDAMB231_trinity.Trinity.fasta.transdecoder.cds -o MDAMB231_trinity.Trinity.fasta.transdecoder.pep

python search_seq_in_assembly_fuzzy.py -i HCC1954_trinity.Trinity.fasta.transdecoder.cds -o HCC1954_trinity.Trinity.fasta.transdecoder.pep 

python search_seq_in_assembly_fuzzy.py -i MDAMB453_trinity.Trinity.fasta.transdecoder.cds -o MDAMB453_trinity.Trinity.fasta.transdecoder.pep 

"""
