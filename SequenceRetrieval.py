#!/usr/bin/env python3

############################################################
# Script for retrieving known stylin multifasta from NCBI  #
############################################################

from Bio import SeqIO
from Bio import Entrez

# All stylins identified and their NCBI ID
# 01, 02, 03, 04, 04bis, 05
NCBI_Gene_ID = ["NM_001162314.1",
                "NM_001162671.2",
                "NM_001161959.2",
                "NM_001172260.1",
                "NM_001172268.2",
                "NM_001163252.1"
                ]
NCBI_Protein_ID = ["NP_001155786.1",
                    "NP_001156143.1",
                    "NP_001155431.1",
                    "NP_001165731.1",
                    "NP_001165739.1",
                    "NP_001156724.1"
                    ]

Entrez.email = "elian.strozyk@etu.umontpellier.fr"

################ Retrieve DNA sequences
dna_search_results = Entrez.efetch("Nucleotide",
        id = ",".join(NCBI_Gene_ID),
        rettype = "fasta")
dna_records = SeqIO.parse(dna_search_results, "fasta")
SeqIO.write(dna_records, "./Stylin_DNA_sequence.fasta", "fasta")

################ Retrieve protein sequences
protein_search_results = Entrez.efetch("Protein",
        id = ",".join(NCBI_Protein_ID),
        rettype = "fasta")
protein_records = SeqIO.parse(protein_search_results, "fasta")
SeqIO.write(protein_records, "./Stylin_protein_sequence.fasta", "fasta")

