#!/usr/bin/env python
# coding: utf-8
from fuzzysearch import find_near_matches
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
import collections
# Load the Pandas libraries with alias 'pd' 
import pandas as pd 
# open file to read
import numpy as np

data = pd.read_csv("cell_lines.txt", sep="\t") #, skiprows = 1)

#print(data.head())
#data_trans = data.T
#print(data_trans)
#np.savetxt(r'np.txt', data_trans.values, fmt='%s')
#data_trans.to_csv(r'pandas.txt',sep='\t', mode='a')
#data_dict = data_trans.to_dict()
#print(data_dict)
cell_dict = defaultdict(list)

f = open('pandas.txt', "r")
for line in f:
   line = line.rstrip()
   if line.split():
      element = line.split("\t")
      #print(element)
      cell_dict[line.split("\t")[0]] += element


#print(cell_dict)

#print(cell_dict["RPMI_8226_EV"])

f_in = open("Tantigen.txt", "r")
line_dict = defaultdict(str)
collect_hits = set([])
count = 0
for line in f_in:
    count = count + 1
    #print(count)
    cloumn1 = line.split()[0]
    if line.startswith("#"):
        continue
    if cloumn1 == "":
        continue
    if len(line.split()) > 1:
        colA = line.split()[0]
        colA = colA.rstrip()

        for cell_line, peptide_list in cell_dict.items():
           for peptide in peptide_list:
               peptide = peptide.rstrip()
               if peptide in colA or colA in peptide:
                   result = "".join([cell_line, "\t", colA, "\t", peptide, "\t"])
                   #if not result in collect_hits:
                   out = "\t".join([cell_line, colA, peptide, "\t"])
                       #out = count, "\t", colA, "\t", temp, " (test peptide)", "\t", data
                   collect_hits.add(result)
                   line_dict[count] += out
for i in range(0, 11403):
    vals = line_dict[i]
    print(i, "\t", vals)
