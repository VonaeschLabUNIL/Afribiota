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
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

cd /scratch/syersin2/Satellite_scratch/output_data/METAPHLAN
mkdir MERGED

merge_metaphlan_tables.py *_metaphlan_profile.txt > ./MERGED/Satellite_metaphlan_relab.txt
merge_metaphlan_tables_abs.py *_metaphlan_profile.txt > ./MERGED/Satellite_metaphlan_abs.txt

cd /scratch/syersin2/Satellite_scratch/output_data/METAPHLAN/GTDB_tax
merge_metaphlan_tables.py --gtdb_profiles *_gtdb.txt > /scratch/syersin2/Satellite_scratch/output_data/METAPHLAN/MERGED/Satellite_metaphlan_gtdb.txt