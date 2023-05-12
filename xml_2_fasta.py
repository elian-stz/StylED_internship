#!/usr/bin/env python3

###########################################################
# Parse an XML file and converts it into a fasta file     #
# Input: an XML file                                      #
# Output: a fasta file                                    #
###########################################################

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys, re 

with open("input.xml") as file :
    for line in file :
        print(line)
