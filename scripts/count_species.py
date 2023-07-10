import sys
from Bio import SeqIO

def counter(fasta_file, type_seq) :
    # create a species_count dictionary: id(key); number of occurrences(value)
    species_count = {}
    with open(fasta_file, "r") as file :
        for record in SeqIO.parse(file, "fasta") :
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
            if species not in species_count :
                species_count[species] = 1
            else :
                species_count[species] += 1
    
    # sorting the species_count by ascending values
    species_count = dict(sorted(species_count.items(), key = lambda item: item[1]))
                
    # printing the species_count
    for key in species_count :
        print(f"{key} : {species_count[key]}")
    
    # Storing the total number of species
    speciesTotalNumber = sum(species_count.values())
    print(speciesTotalNumber)

if __name__ == "__main__" :
    counter(sys.argv[1], sys.argv[2])