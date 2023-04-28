#!/usr/bin/env python3

from Bio.Blast import NCBIWWW
from Bio import SeqIO

with open("Stylin_protein_sequence.fasta", "r") as fd :
    with open("Stylin_protein_blast.xml", "w") as file :
        results = NCBIWWW.qblast("blastp", "nr", fd, 
                entrez_query = 'txid27482[ORGN]' 
        )
        file.write(results.read())
