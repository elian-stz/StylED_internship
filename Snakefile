import os
import pandas as pd
import datetime

#file_list_IDs = "stylin_NCBI_IDs_short.csv" # stylin_NCBI_IDs.csv
#mode = "tblastn"
#organism = "Aphidini" # Aphidini for shorter query
local = False

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

with open(file_list_IDs, "r") as csvf :
	df = pd.read_csv(csvf, sep = ",")
	if mode == "blastp" or mode == "tblastn" :
		ids = df["ncbi_protein_id"].tolist()
	elif mode == "blastn" or mode == "blastx" :
		ids = df["ncbi_gene_id"].tolist()

rule all:
	input:
		expand("{mode}_{organism}/Tree_per_orthogroup", mode=mode, organism=organism),
		expand("{mode}_{organism}/Sequence_number_per_species.txt", mode=mode, organism=organism),
		expand("{mode}_{organism}/Summary_per_orthogroup.csv", mode=mode, organism=organism),
		expand("{mode}_{organism}/Summary_per_orthogroup_selected.csv", mode=mode, organism=organism)

rule run_blast:
	output:
		directory = directory(expand("{mode}_{organism}/Fasta_per_input_ID/", mode=mode, organism=organism)),
		fasta_files = expand("{mode}_{organism}/Fasta_per_input_ID/{id}_{mode}_output.fasta", mode=mode, organism=organism, id=ids)
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
		temp(expand("{mode}_{organism}/Fasta_raw/temp", mode=mode, organism=organism))
	shell:
		"""
		cat {input} > {output}
		"""

rule remove_redundant_sequences:
	input:
		rules.concatenate_fasta_files.output	
	output:
		expand("{mode}_{organism}/Fasta_raw/whole_{mode}_output.fasta", mode=mode, organism=organism)
	shell:
		"""
		python3 scripts/remove_redundant_sequences.py --i {input} --o {output}
		"""

rule split_fasta_files_by_species:
	input:
		rules.remove_redundant_sequences.output
	output:
		directory(expand("{mode}_{organism}/Fasta_per_species/", mode=mode, organism=organism))
	params:
		seq_type = seq_type
	shell:
		"""
		python3 scripts/split_fasta_by_species.py --i {input} --type {params} --outdirectory {output}
		"""

checkpoint run_orthofinder:
	input:
		rules.split_fasta_files_by_species.output
	output:
		orthogroup_dir = directory(expand("{mode}_{organism}/Fasta_per_orthogroup", mode=mode, organism=organism)),
		orthogroups_file = expand("{mode}_{organism}/Orthogroups.txt", mode=mode, organism=organism)
	params:
		seq_type = seq_type,
		orthofinder = "scripts/OrthoFinder/orthofinder",
		orthogroup_dir = directory(expand("{mode}_{organism}/Fasta_per_species/OrthoFinder/Results_{date}/Orthogroup_Sequences", mode=mode, date=current_date, organism=organism)),
		OrthogroupsTxt = expand("{mode}_{organism}/Fasta_per_species/OrthoFinder/Results_{date}/Orthogroups/Orthogroups.txt", mode=mode, organism=organism, date=current_date),
		dir = directory(expand("{mode}_{organism}", mode=mode, organism=organism))
	shell:
		"""
		rm -rf {output.orthogroup_dir}
		if [[ {params.seq_type} = "nucleotide" ]]; then
			{params.orthofinder} -d -f {input}
		elif [[ {params.seq_type} = "protein" ]]; then
			{params.orthofinder} -f {input}
		fi
		cp -r {params.orthogroup_dir} {output.orthogroup_dir}
		cp {params.OrthogroupsTxt} {params.dir}
		"""

rule run_mafft:
	input:
		rules.run_orthofinder.output.orthogroup_dir
	output:
		directory(expand("{mode}_{organism}/Fasta_per_orthogroup_aligned", mode=mode, orthogroup="{orthogroup}", organism=organism))
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
		directory(expand("{mode}_{organism}/IQ-TREE_output_per_orthogroup", mode=mode, organism=organism))	
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
		directory(expand("{mode}_{organism}/Tree_per_orthogroup", mode=mode, organism=organism))
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
		expand("{mode}_{organism}/Sequence_number_per_species.txt", mode=mode, organism=organism)
	params:
		seq_type = seq_type,
		outdir = directory(expand("{mode}_{organism}/", mode=mode, organism=organism))
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
		expand("{mode}_{organism}/Summary_per_orthogroup.csv", mode=mode, organism=organism),
		expand("{mode}_{organism}/Summary_per_orthogroup_selected.csv", mode=mode, organism=organism)
	params:
		seq_type = seq_type,
		csv = file_list_IDs,
		dir = directory(expand("{mode}_{organism}", mode=mode, organism=organism))
	shell:
		"""
		python3 scripts/summarise_orthogroup_as_csv.py --fasta {input.fasta_guide} --type {params.seq_type} --txt {input.orthogroups} --csv {params.csv}
		mv Summary_per_orthogroup.csv {params.dir}
		mv Summary_per_orthogroup_selected.csv {params.dir}
		rm {input.orthogroups}
		"""