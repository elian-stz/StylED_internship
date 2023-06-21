#!/usr/bin/env python3

import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez

def parseInput(infile, mode) :
    # parse a CSV file to retrieve an ID list based on blast algorithm
    with open(infile, "r") as file :
        df = pd.read_csv(file, sep = ",")
        if mode == "blastp" or mode == "tblasn" :
            return df["ncbi_protein_id"].tolist()
        elif mode == "blastn" or mode == "blastx" :
            return df["ncbi_gene_id"].tolist()

def blast(list_id, mode, organism) :
    # Entrez information
    Entrez.email = "elian.strozyk@etu.umontpellier.fr"

    # blast parameters
    match mode :
        case "blastp" :
            database = "nr"
            word_size = 3
            matrix_name = "BLOSUM62"
            gapcosts = "11 1"
            nucl_penalty = None
            nucl_reward = None
            type_returned = "Protein"
        case "tblastn" :
            database = "nr"
            word_size = 3
            matrix_name = "BLOSUM62"
            gapcosts = "11 1"
            nucl_penalty = None
            nucl_reward = None
            type_returned = "Nucleotide"
        case "blastn" :
            database = "nt"
            word_size = 16
            matrix_name = None
            gapcosts = None
            nucl_penalty = -2
            nucl_reward = 1 
            type_returned = "Nucleotide"
        case "blastx" :
            database = "nr"
            word_size = 3
            matrix_name = "BLOSUM62"
            gapcosts = None
            nucl_penalty = None
            nucl_reward = None
            type_returned = "Protein"

    # create a directory to store files
    directory = f"./NCBI_{mode}/Fasta_per_input_ID"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # run blast
    id_seqIO = {} # contains all IDs per input ID 
    for id in list_id :
        result_handle = NCBIWWW.qblast(program = mode, 
                                       database = database,
                                       sequence = id,
                                       entrez_query = organism + "[ORGN]",
                                       expect = 0.5,
                                       gapcosts = gapcosts,
                                       hitlist_size = 500,
                                       word_size = word_size,
                                       matrix_name = matrix_name,
                                       nucl_penalty = nucl_penalty,
                                       nucl_reward = nucl_reward
                                      )
        for record in NCBIXML.parse(result_handle) :
            for alignment in record.alignments :
                if id not in id_seqIO :
                    id_seqIO[id] = []
                id_seqIO[id].append(alignment.hit_id)
        result_handle.close()

    # retrieve sequences as fasta file per input ID
    for key in id_seqIO :
        search_results = Entrez.efetch(type_returned,
                                       id = ",".join(id_seqIO[key]),
                                       rettype = "fasta"
                                      )
        records = SeqIO.parse(search_results, "fasta")
        SeqIO.write(records, f"{directory}/{key}_{mode}_output.fasta", "fasta")
        search_results.close()

def main() :
    if len(sys.argv) != 7 :
        print("error")
        print(f"python3 {sys.argv[0]} --i <csv_file> --algorithm <blast_algorithm_type> --organism <Organism>")
    else :
        list = parseInput(sys.argv[2], sys.argv[4])
        blast(list, sys.argv[4], sys.argv[6])

if __name__ == "__main__" :
    main()
