#!/usr/bin/env python3

###########################################################
# Parses a multifasta file and breaks it down into        #
# several fasta files per species (1 file = 1 species)    #
# This scripts works with fasta files following           #
# NCBI's standard                                         #
#                                                         #
# Input:  a fasta file                                    #
# Output: a directory with several fasta files            #
###########################################################
"""
This script will only work with fasta files following NCBI's standards:
    NCBI header for nucleotide:
        <ID> <species> <name>, <other_info>
        <ID> PREDICTED: <species> <name>, <other_info>
    NCBI header for protein:
        <ID> <name> [<species>] 
"""

from Bio import SeqIO # SeqIO.parse()
import sys            # sys.argv[]
import os             # os.path.exists(), os.makedirs()

def splitter(infile, sequence_type, out_directory) :
    # lowercase for the {sequence_type}
    sequence_type = sequence_type.lower()

    # reuse the name of the input file for the new files
    file_name = infile.split("/")[-1]
    if "whole_" in infile :
        file_name = file_name.replace("whole_", "") # applies only to the pipeline

    # create a directory to store fasta files
    if not os.path.exists(out_directory) :
        os.makedirs(out_directory)

    # preprocess sequences: sort by species
    dictionary = {} # stores SeqIO objects per species
    with open(infile, "r") as file :
        for record in SeqIO.parse(file, "fasta") :
            header = record.description.split() # split after each space
            # NCBI header for nucleotide: <ID> (PREDICTED:) <species> <name>, <other_info>
            # NCBI header for protein:    <ID> <name> [<species>] 
            if sequence_type == "p" or sequence_type == "protein" : 
                species = " ".join(header[-2:])
            if sequence_type == "nt" or sequence_type == "n" or sequence_type == "nucleotide" : 
                if header[1] == "PREDICTED:" or header[1] == "predicted:" :
                    species = " ".join(header[2:4])
                else :
                    species = " ".join(header[1:3])
            # fill the dictionary
            if species not in dictionary :
                dictionary[species] = [record] # list of SeqIO objects
            else :
                dictionary[species].append(record)
        # write the files
        for species, sequences in dictionary.items() :
            # change the species' name to be file-ready 
            species_clean = species.replace(" ", "_").replace("[", "").replace("]", "")
            SeqIO.write(sequences, f"{out_directory}/{species_clean}_{file_name}", "fasta")

def main() :
    if len(sys.argv) != 7 :
        print("This script parses a multifasta file and breaks it down into fasta files per species.")
        print("Each file contains sequences from one species (1 file = 1 species).")
        print("The script assumes the fasta's headers follow NCBI's standard.")
        print("Input: a fasta file")
        print("Output: a directory with several fasta files")
        print("To use this script, type:")
        print(f"python3 {sys.argv[0]} --i <fasta_file> --type <nucleotide/protein> --outdirectory <output_directory>")
    else :
        splitter(sys.argv[2], sys.argv[4], sys.argv[6])

if __name__ == "__main__" :
    main()
