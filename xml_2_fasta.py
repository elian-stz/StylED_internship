#!/usr/bin/env python3

###########################################################
# Parse an XML file and converts it into a fasta file     #
# Input: an XML file, an output file name                 #
# Output: a fasta file (named after the output file name) #
###########################################################

from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys 

def converter(input_file, output_file) : 
    records = []
    with open(input_file, "r") as finput :
        for blast_record in NCBIXML.parse(finput) :
            print("A sequence")
            #for i in blast_record.descriptions :
                #print(i.title)
            for alignment in blast_record.alignments :
                print(alignment.title)
                for hsp in alignment.hsps:
                    seq_record = SeqRecord(Seq(hsp.sbjct), id = alignment.hit_id, description = "")
                    records.append(seq_record)
    with open(output_file, "w") as foutput :
        SeqIO.write(records, foutput, "fasta")

def main() :
    if len(sys.argv) != 5 :
        print("This script parses a BLAST's XML file and converts it into a fasta file.")
        print("Input: an XML file (obtained from a BLAST query)")
        print("Output: a fasta file")
        print("")
        print("To use this script, type:")
        print("./xml_2_fasta.py --i <XML file> --o <fasta file>")
        print("or")
        print("python3 xml_2_fasta.py --i <XML file> --o <fasta file>")
    else :
        converter(sys.argv[2], sys.argv[4])
        print("Conversion succesfully completed.")
        print("The output is located in ", sys.argv[4])
		
if __name__ == "__main__" : 
    main()
