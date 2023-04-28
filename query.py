#!/usr/bin/env python3

from Bio import Entrez
from Bio import SeqIO

Entrez.email = "elian.strozyk@etu.umontpellier.fr"

query = Entrez.esearch(db = "Genome", term = "(Arthropoda[Organism])") 
result = Entrez.read(query)

for i in result :
    print(i, " : ", result[i])

query.close()
list_id = result["IdList"]
sequences = Entrez.efetch(db = "Genome", id = list_id, rettype = "fasta")
sequences.close()
