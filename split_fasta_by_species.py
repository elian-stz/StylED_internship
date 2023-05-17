#!/usr/bin/env python3

###########################################################
# Parse a multifasta file and break it down to several    #
# files with only one species                             #
# Input: a fasta file                                     #
# Output: fasta files                                     #
###########################################################

from Bio import SeqIO
import sys
import re

def splitter(istream) :
    dictionary = {}
    with open(istream, "r") as file :
        for record in SeqIO.parse(istream, "fasta") :
            match = re.search("\[([A-Z][a-z]+ [a-z]+)\]", record.description)
            species = match.group(0)
            if species not in dictionary :
                dictionary[species] = [record]
            else :
                dictionary[species].append(record)
        for species, sequences in dictionary.items() :
            speciesname = species.replace(" ", "_").replace("[", "").replace("]", "")
            SeqIO.write(sequences, f"{speciesname}_{istream}", "fasta")

def main() :
    if len(sys.argv) != 3 :
        print("Parse a multifasta file and breaks it down into several fasta files")
        print("Each fasta contains sequences from the same species")
        print("Input: a fasta file")
        print("Output: several fasta files")
        print("To use this script, type:")
        print("python3 split_fasta_by_species.py --i <fasta_file")
    else :
        splitter(sys.argv[2])
        print("Splitting done.")

if __name__ == "__main__" :
    main()
