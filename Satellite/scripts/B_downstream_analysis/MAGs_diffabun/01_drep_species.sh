#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name drep_species
#SBATCH --output /scratch/syersin2/SATAFB_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/SATAFB_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 50G
#SBATCH --time 12:00:00

# Module
module load miniforge3/4.8.3-4-Linux-x86_64

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/dRep/drep

# Variables
drep_dir=/scratch/syersin2/Satellite_scratch/mags_output/dRep
genomes=/scratch/syersin2/Satellite_scratch/mags_output/MAGs
genome_info=/scratch/syersin2/Satellite_scratch/mags_output/tables

# Dereplicate (species level: 90% completeness and 5% contamination)
dRep dereplicate ${drep_dir} \
    -g ${genomes}/*.fa \
    --genomeInfo ${genome_info}/SAT_checkm_drep_sp.csv \
    -p 32 \
    --S_algorithm fastANI \
    -comp 90 -con 5