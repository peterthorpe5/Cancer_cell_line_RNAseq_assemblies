import sys
import re
import pandas as pd
import numpy as np
from fuzzysearch import find_near_matches
from collections import defaultdict
import collections

##############################################################
# Note this is not my work, but a modified
# verion of https://github.com/lazarlab/XMAn-v2
# updated to newer fasta header format
# this will now use fuzzy matching
##############################################################
# usage: python eliminate_no_mutaa.py protein_id_excel output_file

# obtains mutated peptides from protein description, returns all possible locations of mutated aa in peptide 
def get_peptides(hit):
    peptides = []
    poss_comb = []
    no_hits = hit.split(';')
    for i in range(0, len(no_hits)):
       if no_hits[i]:#.startswith('GN'):
          peptide = no_hits[i].split('|')
          if(len(peptide[3]) == 7):
             peptides.append(peptide[3])
    for pep in range(0,len(peptides)):
       beg = peptides[pep][0:4]
       end = peptides[pep][3:7]
       mid_1 = peptides[pep][1:5]
       mid_2 = peptides[pep][2:6]
       poss_comb = [beg, end, mid_1, mid_2]
    return poss_comb # this is all the combinations
    #return peptides

def main():
   #file names passed as arguments
   excel_file =  sys.argv[1]
   out_file =  sys.argv[2]

   #open and read protein hits excel file; output dataframe
   pd_out = pd.read_excel(excel_file)
   out_df = pd.DataFrame()

   #gets columns sequence and protein description to determine whether mutation is present or not
   sequences = pd_out['Sequence']
   hits = pd_out['Protein Descriptions']
   line_no = 1
   print(pd_out.index)
   for i in range(0, len(pd_out.index)):
      line_no += 1
      if 'GN' in hits[i]:
         peptides = get_peptides(hits[i])
         #print("seq", sequences)
         found = 0
         for x in peptides:
            print(x, "\t", sequences[i])
            MM_s = find_near_matches(x, sequences[i],
                                     max_deletions=1, max_insertions=1,
                                     max_substitutions=2,
                                     max_l_dist=1)
            print(MM_s)
            if MM_s != []:
            #if x in sequences[i]:
               print("yes")
               found = 1
         if found > 0:
            out_df = out_df.append(pd_out.iloc[[i]])
      else:
         out_df = out_df.append(pd_out.iloc[[i]])

   out_df.to_excel(out_file)

main()
