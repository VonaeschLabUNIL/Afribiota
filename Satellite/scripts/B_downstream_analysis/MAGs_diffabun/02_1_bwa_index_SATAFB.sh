#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name index_cat
#SBATCH --output /scratch/syersin2/SATAFB_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/SATAFB_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 50G
#SBATCH --time 08:00:00

# Load necessary modules
module load gcc/12.3.0
module load bwa/0.7.17

# Variables
# Strains
# indir=/scratch/syersin2/SATAFB_scratch/drep/index
# Species
indir=/scratch/syersin2/SATAFB_scratch/drep_species/index

## Index gene catalog
cd ${indir}
# Strains
# bwa index SAT_AFB_all_genomes.fa
# Species
bwa index SAT_AFB_drep_species.fa