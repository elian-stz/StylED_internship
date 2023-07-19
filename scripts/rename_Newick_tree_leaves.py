#!/usr/bin/env python3

###########################################################
# Replace Newick's leaf name (containing only IDs) and    #
# rename them with "species_ID_(commo name)".             #
# Use a guide fasta files containing all IDs in the       #
# Newick file                                             #
# Sequence type and output file name must be typed too    #
#                                                         #
# Inputs: Newick tree                                     #
#         Fasta file, containing IDs in the Newick tree   #
#         Optional: a CSV file, containing common name    #
# Output: Newick tree, with new leaf name                 #
###########################################################

import sys
from Bio import Phylo, SeqIO
import pandas as pd

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
            elif type_seq == "nucleotide" or type_seq == "n" or type_seq == "nt" :
                # NCBI nucleotide fasta header: <ID> (PREDICTED:) <species> <gene_info>
                if "PREDICTED:" in header :
                    species = " ".join(header[2:4])
                else :
                    species = " ".join(header[1:3])
            species = species.replace("[","").replace("]","").replace(" ","_")
            id_nameseq[seq_id] = f"{species}_{seq_id}"
    return id_nameseq

def addCommonNames(dictionary, csv_file, type_seq) :
# append common names contained in a CSV file
    with open(csv_file, "r") as file :
        df = pd.read_csv(file, sep = ",")
        if type_seq == "nucleotide" or type_seq == "nt" or type_seq == "n" :
            ID_common_name = df.set_index("ncbi_gene_id")["common_name"].to_dict()
        if type_seq == "p" or type_seq == "protein" :
            ID_common_name = df.set_index("ncbi_protein_id")["common_name"].to_dict()
    for key in dictionary :
        for id in ID_common_name :
            if key == id :
               dictionary[key] += f"_{ID_common_name[id]}"
    return dictionary

def renameLeaves(newick_file, dictionary, outfile) :
    # rename Newick tree leaves, based on dictionary's key match with ID
    with open(newick_file, "r") as tree_file :
        tree = Phylo.read(tree_file, "newick")
        for leaf in tree.get_terminals() :
            for key in dictionary :
                if key == leaf.name :
                    leaf.name = dictionary[key]
    Phylo.write(tree, outfile, "newick") # write file as Newick
    return 0

def main() :
    if len(sys.argv) == 9 :
        dictionary = processFastaHeaders(sys.argv[4], sys.argv[6])
        renameLeaves(sys.argv[2], dictionary, sys.argv[8])
        return 0
    elif len(sys.argv) == 11 :
        dictionary = processFastaHeaders(sys.argv[4], sys.argv[6])
        dictionary = addCommonNames(dictionary, sys.argv[10], sys.argv[6])
        renameLeaves(sys.argv[2], dictionary, sys.argv[8])
        return 0
    else :
        print("This script renames Newick tree leaves to '{species_name}_{ID_sequence}_({common_name})'")
        print("It assumes that the leaf name is already an ID, contained in a guide fasta file (NCBI format)")
        print("Inputs: Newick file")
        print("        Guide fasta file, containing the IDs in the Newick file")
        print("        Type of sequences: protein or nucleotide")
        print("        Optional: a CSV file containing sequences' common name")
        print("Output: Newick file")
        print("To run this script, type:")
        print(f"python3 {sys.argv[0]} --newick <newick_file> --fasta <guide_fasta_file> --type <nucleotide/protein> --o <output_name>)")
        print("Or, if you have a CSV file:")
        print(f"python3 {sys.argv[0]} --newick <newick_file> --fasta <guide_fasta_file> --type <nucleotide/protein> --o <output_name> --csv <csv_file>")
        return 1

if __name__ == "__main__" :
    main()
