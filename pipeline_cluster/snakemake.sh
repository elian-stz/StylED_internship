#!/bin/bash
#SBATCH -p workq
#SBATCH -t 06:00:00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds". 
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=stefano.colella@inrae.fr
#SBATCH --mail-type=ALL

# Purging modules
module purge

# Loading required modules

module load system/Miniconda3-4.7.10
module load system/Miniconda3
module load bioinfo/OrthoFinder-v2.5.5
module load bioinfo/ncbi-blast-2.6.0+ bioinfo/fastme-2.1.6.1 bioinfo/dlcpar-1.0 bioinfo/mcl-14-137 bioinfo/diamond-v0.9.19
module load bioinfo/mafft-7.487
module load bioinfo/iqtree-2.1.3
module load bioinfo/snakemake-7.8.1
module load system/Python-3.7.4

# Run pipeline 
snakemake -s Snakefile -c2
