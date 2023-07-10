import os
import pandas as pd
import glob

file_list_IDs = "stylin_NCBI_IDs.csv" # stylin_NCBI_IDs.csv		shortIDs.csv
mode = "blastn"
organism = "Aphididae" # Aphidini for shorter query

match mode :
	case "blastp" :
		seq_type = "protein"
	case "blastn" :
		seq_type = "nucleotide"
	case "tblastn" :
		seq_type = "nucleotide"
	case "blastx" :
		seq_type = "protein"

with open(file_list_IDs, "r") as csvf :
	df = pd.read_csv(csvf, sep = ",")
	if mode == "blastp" or mode == "tblasn" :
		ids = df["ncbi_protein_id"].tolist()
	elif mode == "blastn" or mode == "blastx" :
		ids = df["ncbi_gene_id"].tolist()

rule all:
	input:
		"dummy_orthofinder.txt",
		tree_consensual = expand("{mode}_NCBI/Fasta_raw/tree_complete_cleaned.contree", mode=mode),
		tree_file = expand("{mode}_NCBI/Fasta_raw/tree_complete_cleaned.treefile", mode=mode)


rule run_blast:
	output:
		directory = directory(expand("{mode}_NCBI/Fasta_per_input_ID/", mode=mode)),
		fasta_files = expand("{mode}_NCBI/Fasta_per_input_ID/{id}_{mode}_output.fasta", mode=mode, id=ids)
	params:
		list_id = ",".join(ids),
		mode = mode,
		organism = organism
	shell:
		"""
		python3 scripts/blast.py --id {params.list_id} --algorithm {params.mode} --organism {params.organism} --outdirectory {output.directory}
		"""

rule concatenate_fasta_files:
	input:
		files = expand("{mode}_NCBI/Fasta_per_input_ID/{id}_{mode}_output.fasta", mode=mode, id=ids)
	output:
		outfile = expand("{mode}_NCBI/Fasta_raw/temp", mode=mode)
	shell:
		"""
		cat {input} > {output}
		"""

rule remove_redundant_sequences:
	input:
		redundant_file = expand("{mode}_NCBI/Fasta_raw/temp", mode=mode)
	output:
		nonredundant_file = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output.fasta", mode=mode)
	shell:
		"""
		python3 scripts/remove_redundant_sequences.py --i {input} --o {output}
		rm {input}
		"""

rule split_fasta_files_by_species:
	input:
		files = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output.fasta", mode=mode)
	output:
		outdirectory = directory(expand("{mode}_NCBI/Fasta_per_species/", mode=mode))
	params:
		seq_type = seq_type
	shell:
		"""
		python3 scripts/split_fasta_by_species.py --i {input} --type {params} --outdirectory {output}
		"""

rule run_orthofinder:
	input:
		indirectory = expand("{mode}_NCBI/Fasta_per_species/", mode=mode)
	output:
		dummy = "dummy_orthofinder.txt"
		#directory(expand("{mode}_NCBI/Fasta_per_species/OrthoFinder/", mode = mode))
	params:
		mode = mode,
		command = "scripts/OrthoFinder/orthofinder"
	shell:
		"""
		if [[ {params.mode} = "blastn" || {params.mode} = "tblastn" ]]; then
			{params.command} -d -f {input}
		elif [[ {params.mode} = "blastp" || {params.mode} = "blastx" ]]; then
			{params.command} -f {input}
		fi
		touch {output}
		"""

rule run_mafft:
	input:
		input_file = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output.fasta", mode=mode)
	output:
		aligned_seq = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output_aligned.fasta", mode=mode)
	params:
		mode = mode
	shell:
		"""
		mafft --auto --maxiterate 1000 {input} > {output}
		"""

rule run_iqtree:
	input:
		aligned_seq = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output_aligned.fasta", mode=mode)
	output:
		treefile = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output_aligned.fasta.treefile", mode=mode),
		contree = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output_aligned.fasta.contree", mode=mode)
	params:
		mode=mode
	shell:
		"""
		if [[ {params.mode} = "blastn" || {params.mode} = "tblastn" ]]; then
			iqtree2 -s {input} -mset LG,WAG,JTT -bb 1000
		elif [[ {params.mode} = "blastp" || {params.mode} = "blastx" ]]; then
			iqtree2 -s {input} -m GTR -bb 1000
		fi
		"""

rule change_leaf_name:
	input:
		treefile = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output_aligned.fasta.treefile", mode=mode),
		contree = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output_aligned.fasta.contree", mode=mode),
		fasta = expand("{mode}_NCBI/Fasta_raw/whole_{mode}_output.fasta", mode=mode)
	output:
		tree_consensual = expand("{mode}_NCBI/Fasta_raw/tree_complete_cleaned.contree", mode=mode),
		tree_file = expand("{mode}_NCBI/Fasta_raw/tree_complete_cleaned.treefile", mode=mode)
	params:
		seq_type = seq_type
	shell:
		"""
		python3 scripts/change_leaf_names.py --newick {input.treefile} --fasta {input.fasta} --type {params} --o {output.tree_file}
		python3 scripts/change_leaf_names.py --newick {input.contree} --fasta {input.fasta} --type {params} --o {output.tree_consensual}
		"""
