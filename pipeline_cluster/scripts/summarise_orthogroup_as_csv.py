#!/usr/bin/env python3

###########################################################
# Parse OrthoFinder's output file named Orthogroups.txt   #
# obtained after running OrthoFinder. It summarises data  #
# from fasta's header as a csv file for biologists.		  #
#														  #
#  Input: a txt file (named Orthogroups.txt)			  #
# Output: two csv files									  #
#				one containing all orthogroups			  #
#				one containing selected orthogroups		  #
###########################################################

from Bio import SeqIO
import pandas as pd
import sys 

def processFastaHeaders(fasta_file, seq_type) :
	# create a dictionary: id(key); fasta header split (value)
	id_header = {}
	with open(fasta_file, "r") as file :
		for record in SeqIO.parse(file, "fasta") :
			seq_id = record.id
			header = record.description.split()
			if seq_type == "protein" or seq_type == "p" :
				# NCBI protein fasta header: <ID> <protein_name> [<species>]
				species = " ".join(header[-2:])
				name = " ".join(header[1:-2])
			elif seq_type == "nucleotide" or seq_type == "n" or seq_type == "nt" :
				# NCBI nucleotide fasta header: <ID> (PREDICTED:) <species> <gene_info>
				if "PREDICTED:" in header :
					species = " ".join(header[2:4])
					name = "PREDICTED: " + " ".join(header[4:])
				else :
					species = " ".join(header[1:3])
					name = " ".join(header[3:])
			species = species.replace("[","").replace("]","")
			id_header[seq_id] = [seq_id, species, name]
	return id_header

def parseOrthogroupsFile(txtfile, dictionary) :
# modify the dictionary previously created by appending the orthogroup to the value
	with open(txtfile, "r") as file :
		for line in file :
			element = line.split() # 1st element is orthogroup, others are IDs
			orthogroup = element[0].replace(":","")
			for seq_id in element[1:] :
				for key in dictionary :
					if key == seq_id :
						dictionary[key].append(orthogroup)
	return dictionary

def retrieveInputIDs(csvfile, seq_type) :
# retrieve from the CSV file the corresponding IDs depending on the sequence type
	with open(csvfile, "r") as file :
		df_guide = pd.read_csv(file, sep = ",")
		if seq_type == "nucleotide" or seq_type == "nt" or seq_type == "n" :
			input_IDs = df_guide["ncbi_gene_id"].tolist()
		if seq_type == "p" or seq_type == "protein" :
			input_IDs = df_guide["ncbi_protein_id"].tolist()
	return input_IDs

def makeOutput(dictionary, list) :
	#convert dictionary's values into list of lists and add a dataframe header
	sequence_headers = []
	for key in dictionary :
		sequence_headers.append(dictionary[key])
	csv_header = ["NCBI_ID", "Species", "Sequence_name", "Orthogroup"]
	sequence_headers.insert(0, csv_header)

	# convert into dataframe
	df_main = pd.DataFrame(sequence_headers[1:], columns = sequence_headers[0])

	# add an attribute: if NCBI_ID row contains the identified IDs, add a 1(True)
	df_main["Identified_input_ID"] = df_main["NCBI_ID"].isin(list).astype(int)

	# reorder attributes
	df_main = df_main.reindex(columns = ["Orthogroup",
										 "Species",
										 "Identified_input_ID", 
										 "NCBI_ID",
										 "Sequence_name"
										]
							 )

	# write the file
	df_main.to_csv("Summary_per_orthogroup.csv", index = False)

	#create a subfile containing only orthogroups with identified input ID 
	groups = df_main.loc[df_main["Identified_input_ID"] == 1, "Orthogroup"]
	df_selected = df_main[df_main["Orthogroup"].isin(groups.tolist())]
	df_selected.to_csv("Summary_per_orthogroup_selected.csv", index = False) # write the subfile 

def main() :
	if len(sys.argv) != 9 :
		print("This script parses OrthoFinder's output file named Orthogroups.txt")
		print("It returns two csv files summarising header's fasta")
		print("The input file is located in /Results_<date>/Orthogroups/ after running OrthoFinder")
		print("Input: a txt file (named Orthogroups.txt)")
		print("Output: two csv files")
		print("		   one containing all orthogroups")
		print("		   one containing only selected orthogroups with proteins of interest (identified stylins)")
		print("To run this script, type:")
		print(f"python3 {sys.argv[0]} --fasta <guide_fasta> --type <nucleotide/protein> --txt <Orthogroups.txt> --csv <guide_csv>")
	else :
		dictionary = processFastaHeaders(sys.argv[2], sys.argv[4])
		dictionary = parseOrthogroupsFile(sys.argv[6], dictionary)
		inputIDs = retrieveInputIDs(sys.argv[8], sys.argv[4])
		makeOutput(dictionary, inputIDs)

if __name__ == "__main__" :
	main()
