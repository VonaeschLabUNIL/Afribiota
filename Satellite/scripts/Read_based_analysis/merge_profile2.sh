#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name merge_metaphlan
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 8G
#SBATCH --time 00:15:00
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err

# Module
module load gcc/13.2.0
module load miniforge3/4.8.3-4-Linux-x86_64

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

cd /scratch/syersin2/metaphlan/profiles
mkdir MERGED

merge_metaphlan_tables.py *.txt > ./MERGED/AFB_SAT_metaphlan_relab.txt
merge_metaphlan_tables_abs.py *.txt > ./MERGED/AFB_SAT_metaphlan_abs.txt