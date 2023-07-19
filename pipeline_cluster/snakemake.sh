#!/bin/bash
#SBATCH -p workq
#SBATCH -t 06:00:00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds". 
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=elian.strozyk@etu.umontpellier.fr
#SBATCH --mail-type=ALL

# Loading required modules
module load bioinfo/snakemake-7.8.1
module load system/Python-3.7.4
module load bioinfo/ncbi-blast-2.10.0+ 
module load bioinfo/fastme-2.1.6.1 
module load bioinfo/dlcpar-1.0
module load bioinfo/mcl-14-137 
module load bioinfo/diamond-v0.9.19
module load bioinfo/OrthoFinder-2.5.4
module load bioinfo/mafft-7.487
module load bioinfo/iqtree-2.1.3

# Run pipeline 
snakemake -s Snakefile --cores 1 --keep-going --cluster-config SCRIPTS/cluster_config.yaml --cluster "sbatch -J {cluster.job-name} -A recombinationlandscape -p {cluster.queue} {cluster.nodes} {cluster.cpus} {cluster.log} {cluster.mem}" --jobs 5
