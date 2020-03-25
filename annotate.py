from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
from collections import defaultdict
import os
from sys import stdin,argv


def seq_getter(transcript_pep, blast, outfile):
    """fun to add the desc from UP000005640_9606.fasta to the
    top blast hits"""
    f= open(outfile, 'w')

    human_prot =  SeqIO.index("UP000005640_9606.fasta", "fasta")

    blast_open = open(blast, "r")
    trans_to_human = defaultdict(str)
    for line in blast_open:
        if line.strip():
            line = line.rstrip()
            data = line.split("\t")
            transcript = data[0]
            human_hit = data[1]
            human_record = human_prot[human_hit]
            outfmt = "%s" % (human_record.description)
            trans_to_human[transcript] = outfmt
    blast_open.close()
    for seq_record in SeqIO.parse(transcript_pep, "fasta"):
        new_desc = trans_to_human[seq_record.id]
        seq_record.description = new_desc
        SeqIO.write(seq_record, f, "fasta")
  
    f.close()


seq_getter(argv[1],argv[2], argv[3])
