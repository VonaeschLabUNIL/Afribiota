#!/bin/bash

#SBATCH --mem 10G 
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 08:00:00
#SBATCH --cpus-per-task 24
#SBATCH --job-name phylogeny
#SBATCH --output /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.err

# Module
module load gcc
module load miniconda3

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/OrthoFinder/orthofinder

# Variables
workDir=/scratch/syersin2/Satellite_scratch/Isolates/phylogeny/OrthoFinder

## CHANGE DATE!!
data=${workDir}/Results_Dec11/MultipleSequenceAlignments

cd ${workDir}

iqtree2 -s ${data}/SpeciesTreeAlignment.fa\
        -st AA -nt 16 -bb 1000 -seed 12345 -m TEST\
        -pre SAT_isolates
