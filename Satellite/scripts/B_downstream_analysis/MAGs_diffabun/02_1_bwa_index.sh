#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name index_cat
#SBATCH --output /scratch/syersin2/Satellite_scratch/mags_output/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/mags_output/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 50G
#SBATCH --time 08:00:00

# Load necessary modules
module load gcc/12.3.0
module load bwa/0.7.17

# Variables
# Species
indir=/scratch/syersin2/Satellite_scratch/mags_output/dRep/index

## Index gene catalog
cd ${indir}
# Species
bwa index Satellite_drep_species.fa