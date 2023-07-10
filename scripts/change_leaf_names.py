#!/usr/bin/env python3

###########################################################
# Replace Newick's leaf names (containing only IDs) and   #
# rename them with "species_ID". Use a guide fasta files  #
# containing all IDs in the Newick file                   #
# Sequence type and output file name must be typed too    #
#                                                         #
# Inputs: Newick tree                                     #
#         Fasta file, containing IDs in the Newick tree   #
# Output: Newick tree, with new leaf name                 #
###########################################################

import sys
from Bio import Phylo, SeqIO
import re

# newick_file: tree to change name
# fasta file to retrieve header
# output_file : new file

def processFastaHeaders(fasta_file, type_seq) :
    # create a dictionary: id(key); new leaf name (value)
    id_nameseq = {}
    with open(fasta_file, "r") as file :
        for record in SeqIO.parse(file, "fasta") :
            seq_id = record.id
            header = record.description.split()
            if type_seq == "protein" or type_seq == "p" :
                # NCBI protein fasta header: <ID> <protein_name> [<species>]
                species = " ".join(header[-2:])
            elif type_seq == "nucleotide" or type_seq == "n" :
                # NCBI nucleotide fasta header: <ID> (PREDICTED:) <species> <gene_info>
                if "PREDICTED:" in header :
                    species = " ".join(header[2:4])
                else :
                    species = " ".join(header[1:3])
            species = species.replace("[","").replace("]","").replace(" ","_")
            id_nameseq[seq_id] = f"{species}_{seq_id}"
    return id_nameseq

def changeLeafName(newick_file, dictionary, outfile) :
    # replace leaf name using a dictionary, based on dictionary's key match with ID
    with open(newick_file, "r") as tree_file :
        tree = Phylo.read(tree_file, "newick")
        for leaf in tree.get_terminals() :
            for key in dictionary :
                if key == leaf.name :
                    leaf.name = dictionary[key]
    Phylo.write(tree, outfile, "newick") # write file as Newick

def main() :
    if len(sys.argv) != 9 :
        print("This script changes the leaf name of Newick files to '{species_name}_{ID_sequence}'")
        print("It assumes that the leaf name is already an ID, contained in a guide fasta file (NCBI format)")
        print("Inputs: Newick file")
        print("        Guide fasta file, containing the IDs in the Newick file")
        print("        Type of sequences: protein or nucleotide")
        print("Output: Newick file")
        print("To run this script, type:")
        print(f"python3 {sys.argv[0]} --newick <newick_file> --fasta <guide_fasta_file> --type <nucleotide/protein> --o <output_name>)")
        return 1
    else :
        dictionary = processFastaHeaders(sys.argv[4], sys.argv[6])
        changeLeafName(sys.argv[2], dictionary, sys.argv[8])
        return 0

if __name__ == "__main__" :
    main()
