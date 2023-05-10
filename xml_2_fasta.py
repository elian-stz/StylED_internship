from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys 

def converter(input_file, output_file) : 
	records = []
	with open(input_file, "r") as finput :
		blast_records = NCBIXML.parse(finput)
		for blast_record in blast_records:
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					seq_record = SeqRecord(Seq(hsp.sbjct), id = alignment.hit_id, description = "")
					records.append(seq_record)

	with open(output_file, "w") as foutput :
		SeqIO.write(records, foutput, "fasta")

def main() :
	if len(sys.argv) == 1 :
		print("This script parses a BLAST's XML file and converts it into a fasta file.")
		print("Input: a XML file, obtained from a BLAST query")
		print("Output: a fasta file")
		print("")
		print("To use this script:")
		print("./xml_2_fasta.py --i <XML file> --o <fasta file>")
		print("e.g.: ./xml_2_fasta.py --i blast_results.xml --o blast_sequences.fasta")
	elif len(sys.argv) == 5 :
		converter(sys.argv[2], sys.argv[4])
		
if__name__== "__main__" :
	main()
