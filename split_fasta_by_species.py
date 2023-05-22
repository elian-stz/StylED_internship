#!/usr/bin/env python3

###########################################################
# Parse a multifasta file and breaks it down into several #
# fasta files containing sequences from one species       #
# Input: a fasta file                                     #
# Output: fasta files                                     #
###########################################################

from Bio import SeqIO
import sys
import re

def splitter(infile) :
    dictionary = {} # stores SeqIO objects per species
    with open(infile, "r") as file :
        for record in SeqIO.parse(file, "fasta") :
            # assuming the fasta's header contains the species written like this: [Myzus persicae]
            match = re.search("\[([A-Z][a-z]+ [a-z]+)\]", record.description)
            species = match.group(0)
            if species not in dictionary :
                dictionary[species] = [record] # list of SeqIO objects
            else :
                dictionary[species].append(record)
        for species, sequences in dictionary.items() :
            # change the species' name to be file-ready 
            species_file = species.replace(" ", "_").replace("[", "").replace("]", "")
            SeqIO.write(sequences, f"{species_file}_{file}", "fasta") # write the file

def main() :
    if len(sys.argv) != 3 :
        print("Parse a multifasta file and breaks it down into several fasta files")
        print("Each file contains sequences from one species")
        print("The script assumes the fasta's headers contains this pattern: [Homo sapiens]")
        print("Input: a fasta file")
        print("Output: several fasta files")
        print("To use this script, type:")
        print("python3 split_fasta_by_species.py --i <fasta_file>")
    else :
        splitter(sys.argv[2])
        print("Splitting done.")

if __name__ == "__main__" :
    main()
