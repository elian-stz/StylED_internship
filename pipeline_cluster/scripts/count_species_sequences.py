import sys
from Bio import SeqIO

def countSpecies(fasta_file, type_seq) :
	# create a species_count dictionary: species(key); number of occurrences(value)
	species_count = {}
	with open(fasta_file, "r") as file :
		for record in SeqIO.parse(file, "fasta") :
			header = record.description.split()
			if type_seq == "protein" or type_seq == "p" :
			# NCBI protein fasta header: <ID> <protein_name> [<species>]
				species = " ".join(header[-2:])
			elif type_seq == "nucleotide" or type_seq == "n" :
			# NCBI nucleotide fasta header: <ID> (PREDICTED:) <species> <gene_name>, <gene_info>
				if "PREDICTED:" in header :
					species = " ".join(header[2:4])
				else :
					species = " ".join(header[1:3])
			if species not in species_count :
				species_count[species] = 1
			else :
				species_count[species] += 1
	return species_count

def writeOutput(dictionary) :
	with open("Sequence_number_per_species.txt", "w") as file:
		file.write(f"Total number of sequences: {sum(dictionary.values())}\n")
		for key in dictionary :
			file.write(f"{key}: {dictionary[key]}\n")

def main() :
	if len(sys.argv) != 5 :
		print("This script generates a file that sums up the number of sequences per species.")
		print("It assumes that fasta headers follow NCBI's standards.")
		print("Input: a fasta file")
		print("Output: a fasta file")
		print("To run this script, type:")
		print(f"python3 {sys.argv[0]} --i <fasta_file> --type <nucleotide/protein>")
	else:
		dictionary = countSpecies(sys.argv[2], sys.argv[4])
		writeOutput(dictionary)

if __name__ == "__main__" :
	main()
