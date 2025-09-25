#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name SP_strain_share
#SBATCH --output /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 2G
#SBATCH --time 00:05:00
#SBATCH --array=1-480

# Module
module load gcc/12.3.0
module load miniforge3/4.8.3-4-Linux-x86_64

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
indir=/scratch/syersin2/Satellite_scratch/strainphlan/StrainPhlAn_out
outdir=/scratch/syersin2/Satellite_scratch/strainphlan/Strain_sharing
mddir=/users/syersin2/Satellite/scripts/Read_based_analysis

## Define array variable
# Extract the filename corresponding to this SLURM array task ID
cd ${indir}
file=$(ls RAxML_bestTree.t__SGB*.StrainPhlAn4.tre | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract the unique identifier from the filename
SGB=$(basename "${file}" | sed 's/^RAxML_bestTree\.//' | sed 's/\.StrainPhlAn4\.tre$//')

mkdir -p ${outdir}/${SGB}

python /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan/bin/strain_transmission.py \
        --tree ${indir}/RAxML_bestTree.${SGB}.StrainPhlAn4.tre \
        --metadata ${mddir}/cleaned_metadata.tsv \
        --output_dir ${outdir}/${SGB} 