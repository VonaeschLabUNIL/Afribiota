#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name CheckM
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 36
#SBATCH --mem 200G
#SBATCH --time 02:00:00

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/checkM/checkm

# Variables
indir=/scratch/syersin2/Satellite_scratch/Isolates/wgs_assemblies
outdir=/scratch/syersin2/Satellite_scratch/Isolates/checkm

# Run checkm
mkdir -p ${outdir}

checkm lineage_wf ${indir} \
    ${outdir} \
    -x fasta \
    -t 36 \
    -f ${outdir}/Assemblies_checkM_stats.tsv \
    --tab_table