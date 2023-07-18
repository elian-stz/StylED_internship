configfile: "config.yaml"

import os
import pandas as pd
import datetime

file_list_IDs = config["csv_file"]
mode = config["mode"]
organism = config["organism"]
DIR = f"{mode}_{organism}"

current_date = datetime.date.today().strftime("%b%d")

match mode :
	case "blastp" :
		seq_type = "protein"
	case "blastn" :
		seq_type = "nucleotide"
	case "tblastn" :
		seq_type = "nucleotide"
	case "blastx" :
		seq_type = "protein"

with open(file_list_IDs, "r") as file :
	df = pd.read_csv(file, sep = ",")
	if mode == "blastp" or mode == "tblastn" :
		ids = df["ncbi_protein_id"].tolist()
	elif mode == "blastn" or mode == "blastx" :
		ids = df["ncbi_gene_id"].tolist()

rule all:
	input:
		expand("{DIR}/Tree_per_orthogroup", DIR=DIR),
		expand("{DIR}/Sequence_number_per_species.txt", DIR=DIR),
		expand("{DIR}/Summary_per_orthogroup.csv", DIR=DIR),
		expand("{DIR}/Summary_per_orthogroup_selected.csv", DIR=DIR)

rule run_blast:
	output:
		directory = directory(expand("{DIR}/Fasta_per_input_ID/", DIR=DIR)),
		fasta_files = expand("{DIR}/Fasta_per_input_ID/{id}_{mode}_output.fasta", DIR=DIR, mode=mode, id=ids)
	params:
		list_id = ",".join(ids),
		mode = mode,
		organism = organism,
		local = local
	shell:
		"""
		python3 scripts/blast.py --id {params.list_id} --algorithm {params.mode} --organism {params.organism} --outdirectory {output.directory}
		"""

rule concatenate_fasta_files:
	input:
		rules.run_blast.output.fasta_files	
	output:
		temp(expand("{DIR}/Fasta_raw/temp", DIR=DIR))
	shell:
		"""
		cat {input} > {output}
		"""

rule remove_redundant_sequences:
	input:
		rules.concatenate_fasta_files.output	
	output:
		expand("{DIR}/Fasta_raw/whole_{mode}_output.fasta", DIR=DIR, mode=mode)
	shell:
		"""
		python3 scripts/remove_redundant_sequences.py --i {input} --o {output}
		"""

rule split_fasta_files_by_species:
	input:
		rules.remove_redundant_sequences.output
	output:
		directory(expand("{DIR}/Fasta_per_species/", DIR=DIR))
	params:
		seq_type = seq_type
	shell:
		"""
		python3 scripts/split_fasta_by_species.py --i {input} --type {params} --outdirectory {output}
		"""

rule run_orthofinder:
	input:
		rules.split_fasta_files_by_species.output
	output:
		orthogroup_dir = directory(expand("{DIR}/Fasta_per_orthogroup", DIR=DIR)),
		orthogroups_file = expand("{DIR}/Orthogroups.txt", DIR=DIR)
	params:
		seq_type = seq_type,
		orthofinder = "scripts/OrthoFinder/orthofinder.py",
		orthogroup_dir = directory(expand("{DIR}/Fasta_per_species/OrthoFinder/Results_{date}/Orthogroup_Sequences", DIR=DIR, date=current_date)),
		OrthogroupsTxt = expand("{DIR}/Fasta_per_species/OrthoFinder/Results_{date}/Orthogroups/Orthogroups.txt", DIR=DIR, date=current_date),
		dir = directory(expand("{DIR}", DIR=DIR))
	shell:
		"""
		rm -rf {output.orthogroup_dir}
		if [[ {params.seq_type} = "nucleotide" ]]; then
			python3 {params.orthofinder} -og -d -f {input}
		elif [[ {params.seq_type} = "protein" ]]; then
			python3 {params.orthofinder} -og -f {input}
		fi
		cp -r {params.orthogroup_dir} {output.orthogroup_dir}
		cp {params.OrthogroupsTxt} {params.dir}
		"""

rule run_mafft:
	input:
		rules.run_orthofinder.output.orthogroup_dir
	output:
		directory(expand("{DIR}/Fasta_per_orthogroup_aligned", DIR=DIR))
	shell:
		"""
		mkdir -p {output}
		FILES={input}/*.fa
		for f in ${{FILES}}; do
			if [ -f "$f" ]; then
				mafft --auto --maxiterate 1000 "$f" > "${{f}}.aligned"
			fi
		done
		mv {input}/*.fa.aligned {output}
		"""

rule run_iqtree:
# IQ-TREE doesn't work if the file contains less than 4 sequences
	input:
		rules.run_mafft.output
	output:
		directory(expand("{DIR}/IQ-TREE_output_per_orthogroup", DIR=DIR))	
	params:
		seq_type = seq_type
	shell:
		"""
		mkdir -p {output}
		FILES={input}/*.fa.aligned
		for f in ${{FILES}}; do
			count=$(grep -c ">" "$f")
			if [[ $count -ge 4 ]]; then
				if [[ {params.seq_type} = "protein" ]]; then
					iqtree2 -s "$f" -mset WAG,LG,JTT -bb 1000
				elif [[ {params.seq_type} = "nucleotide" ]]; then
					iqtree2 -s "$f" -m GTR -bb 1000
				fi
			fi
		done
		shopt -s extglob
		mv {input}/!(*.fa.aligned) {output}
		"""

rule rename_tree_leaves:
	input:
		newick_dir = rules.run_iqtree.output,
		fasta = rules.remove_redundant_sequences.output
	output:
		directory(expand("{DIR}/Tree_per_orthogroup", DIR=DIR))
	params:
		seq_type = seq_type,
		csv = file_list_IDs 
	shell:
		"""
		mkdir -p {output}
		CONTREES={input.newick_dir}/*.contree
		TREEFILES={input.newick_dir}/*.treefile
		for f in ${{TREEFILES}} ${{CONTREES}}; do
			python3 scripts/rename_Newick_tree_leaves.py --newick "$f" --fasta {input.fasta} --type {params.seq_type} --o "${{f}}.cleaned" --csv {params.csv} 
		done
		mv {input.newick_dir}/*.cleaned {output}
		"""

rule count_sequences_per_species:
	input:
		rules.remove_redundant_sequences.output
	output:
		expand("{DIR}/Sequence_number_per_species.txt", DIR=DIR)
	params:
		seq_type = seq_type,
		outdir = directory(expand("{DIR}/", DIR=DIR))
	shell:
		"""
		python3 scripts/count_species_sequences.py --i {input} --type {params.seq_type}
		mv Sequence_number_per_species.txt {params.outdir} 
		"""

rule summarise_fasta_files:
	input:
		orthogroups = rules.run_orthofinder.output.orthogroups_file,
		fasta_guide = rules.remove_redundant_sequences.output
	output:
		expand("{DIR}/Summary_per_orthogroup.csv", DIR=DIR),
		expand("{DIR}/Summary_per_orthogroup_selected.csv", DIR=DIR)
	params:
		seq_type = seq_type,
		csv = file_list_IDs,
		dir = directory(expand("{DIR}", DIR=DIR))
	shell:
		"""
		python3 scripts/summarise_orthogroup_as_csv.py --fasta {input.fasta_guide} --type {params.seq_type} --txt {input.orthogroups} --csv {params.csv}
		mv Summary_per_orthogroup.csv {params.dir}
		mv Summary_per_orthogroup_selected.csv {params.dir}
		rm {input.orthogroups}
		"""
