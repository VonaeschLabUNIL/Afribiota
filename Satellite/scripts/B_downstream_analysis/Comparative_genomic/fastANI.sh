#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name fastANI
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G
#SBATCH --time 00:30:00

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/FastANI/fastani

# Variables
genomes=/scratch/syersin2/Satellite_scratch/transmission_dir/genomes
refdir=/scratch/syersin2/Satellite_scratch/comp_genomics/fastANI
outdir=/scratch/syersin2/Satellite_scratch/comp_genomics/fastANI/fastANI_out

mkdir -p ${outdir}

cd ${genomes}
# Run fastANI
fastANI --ql ${refdir}/fastANI_list.txt \
    --rl ${refdir}/fastANI_list.txt \
    -o ${outdir}/SAT_fastANI.tsv \
    -t 8
