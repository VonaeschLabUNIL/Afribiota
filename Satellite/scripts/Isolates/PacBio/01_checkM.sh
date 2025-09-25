#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name CheckM_ISO
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 200G
#SBATCH --time 04:00:00

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/checkM/checkm

# Variables
inputfolder=/scratch/syersin2/Satellite_scratch/Isolates
#outdir=/scratch/syersin2/Satellite_scratch/Isolates/outputs

rm -r ${inputfolder}/outputs/CheckM
mkdir -p ${inputfolder}/outputs/CheckM

# Run checkm
#export CHECKM_DATA_PATH=/path/to/my_checkm_data
checkm lineage_wf ${inputfolder}/assemblies \
    ${inputfolder}/outputs/CheckM \
    -x fasta \
    -t 32 \
    -f ${inputfolder}/outputs/CheckM/Checkm_QC_stats.tsv \
    --tab_table
