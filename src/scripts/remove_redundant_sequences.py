#!/usr/bin/env python3

###########################################################
# Parse a multifasta file and remove redundant sequences  #
# Input: a fasta file                                     #
# Output: a fasta file                                    #
###########################################################

from Bio import SeqIO
import sys

def cleaner(infile, outfile) :
    seen = set() # contains all SeqIO.seq seen
    records = [] # contains all SeqIO objects whose sequences are not in {seen}
    with open(infile, "r") as file :
        for record in SeqIO.parse(file, "fasta") :
            if record.seq not in seen :
                seen.add(record.seq)
                records.append(record)
    SeqIO.write(records, outfile, "fasta") # write file

def main() :
    if len(sys.argv) != 5 :
        print("Parse a multifasta file and remove redundant sequences")
        print("Input: a fasta file")
        print("Output: a fasta file")
        print("To use this script, type:")
        print(f"python3 {sys.argv[0]} --i <fasta_file_to_clean> --o <fasta_file_cleaned>")
    else :
        cleaner(sys.argv[2], sys.argv[4])

if __name__ == "__main__" :
    main()
