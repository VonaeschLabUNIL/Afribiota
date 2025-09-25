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

# Module
module load gcc/12.3.0
module load miniforge3/4.8.3-4-Linux-x86_64

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
indir=/scratch/syersin2/Satellite_scratch/strainphlan/StrainPhlAn_out
outdir=/scratch/syersin2/Satellite_scratch/strainphlan/Strain_sharing_test
mddir=/users/syersin2/Satellite/scripts/Read_based_analysis

## Define array variable
# Extract the filename corresponding to this SLURM array task ID
cd ${indir}

mkdir -p ${outdir}/t__SGB8071

python /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan/bin/strain_transmission.py \
        --tree ${indir}/RAxML_bestTree.t__SGB8071.StrainPhlAn4.tre \
        --metadata ${mddir}/metadata2.tsv \
        --output_dir ${outdir}/t__SGB8071 \
        --threshold 0.03607692237