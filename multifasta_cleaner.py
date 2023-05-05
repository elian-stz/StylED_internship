#!/usr/bin/env python3

###########################################################
# Parse a multifasta file and trim sequences not          #
# belonging to a clade as well as removing redundant      #
# ones                                                    #
#                                                         # 
# INPUT: A SINGLE FASTA FILE                              #
# OUTPUT: A SINGLE FASTA FILE                             #
###########################################################

from Bio import SeqIO
from ete3 import PhyloTree
import re

dictionary = {}
cleaned_seq = []
with open("../metaphors/stylin_orthologs.metaphors.fasta", "r") as fd :
#    for record in SeqIO.parse(fd, format = "fasta") :
#        if record.seq not in cleaned_seq :
#            cleaned_seq.append(record.seq)
#    print(len(cleaned_seq))
	for line in fd :
	        for key in re.findall("\[(.*)\]", line) :
	            if key not in dictionary :
	                dictionary[key] = 1
	            else :
	                dictionary[key] += 1
	for key in dictionary :
        	print(f"{key} : {dictionary[key]}")

