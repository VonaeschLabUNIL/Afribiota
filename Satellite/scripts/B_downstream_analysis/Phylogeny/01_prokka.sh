#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name prokka
#SBATCH --output /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G
#SBATCH --time 00:15:00
#SBATCH --array=1-104

# Module
module load gcc
module load prokka
module load blast-plus

# Prokka is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files.

## Variables
prokka_indir=/scratch/syersin2/Satellite_scratch/Isolates/wgs_assemblies
prokka_outdir=/scratch/syersin2/Satellite_scratch/Isolates/prokka

cd ${prokka_indir}
iso_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p | sed 's/.fasta$//')


echo "#####################################################"
echo "Prokka on isolate": ${iso_name}
echo "#####################################################"

## Run prokka
prokka --outdir ${prokka_outdir}/${iso_name} \
    --prefix ${iso_name} ${prokka_indir}/${iso_name}.fasta \
    --dbdir /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/Prokka/db