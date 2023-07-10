#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez

def blast(list_id, mode, organism, outdirectory) :
    list_id = list_id.split(",")

    # output directory that will store fasta files
    if not os.path.exists(outdirectory):
        os.makedirs(outdirectory)

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
            type_returned = "protein"
        case "tblastn" :
            database = "nr"
            word_size = 3
            matrix_name = "BLOSUM62"
            gapcosts = "11 1"
            nucl_penalty = None
            nucl_reward = None
            type_returned = "nucleotide"
        case "blastn" :
            database = "nt"
            word_size = 16
            matrix_name = None
            gapcosts = None
            nucl_penalty = -2
            nucl_reward = 1 
            type_returned = "nucleotide"
        case "blastx" :
            database = "nr"
            word_size = 3
            matrix_name = "BLOSUM62"
            gapcosts = None
            nucl_penalty = None
            nucl_reward = None
            type_returned = "protein"

    # run blast
    hit_ids = [] # contains all ids retrieved from blast
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
                if alignment.hit_id not in hit_ids:
                    hit_ids.append(alignment.hit_id)
        print(hit_ids)
        result_handle.close()

        # normalise IDs if dual ID e.g. gi|240849398|ref|NM_001162314.1| + remove assembly IDs (OX or OU or OY)
        if type_returned == "nucleotide" :
            cleaned_hit_ids = []
            for hit in hit_ids :
                if "OX" not in hit :
                    if "OU" not in hit :
                        if "OY" not in hit :
                            if len(hit.split('|')) == 5 :
                                cleaned_hit_ids.append(hit.split('|')[3])
                            else :
                                cleaned_hit_ids.append(hit)
            hit_ids = cleaned_hit_ids
        
        # retrieve sequences as fasta file
        search_results = Entrez.efetch(db = type_returned,
                                       id = ",".join(hit_ids),
                                       rettype = "fasta"
                                      )
        records = SeqIO.parse(search_results, "fasta")
        SeqIO.write(records, f"{outdirectory}/{id}_{mode}_output.fasta", "fasta")
        search_results.close()

def main() :
    if len(sys.argv) != 9 :
        print("error")
        print(f"python3 {sys.argv[0]} --id <NCBI_ID_list> --algorithm <blast_algorithm_type> --organism <Organism> --outdirectory <output_directory>")
    else :
        blast(sys.argv[2], sys.argv[4], sys.argv[6], sys.argv[8])

if __name__ == "__main__" :
    main()
