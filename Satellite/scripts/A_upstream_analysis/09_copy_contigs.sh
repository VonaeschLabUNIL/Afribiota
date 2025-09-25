#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name copy_contigs
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 8G
#SBATCH --time 12:00:00
#SBATCH --array=1-129

## Variables
files=/scratch/syersin2/Satellite_scratch/output_data/metaspades

## Array variables
cd ${files}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p /scratch/syersin2/Satellite_scratch/output_data/contigs
cp ${files}/${sample_name}/${sample_name}.contigs.min0.fasta.gz /scratch/syersin2/Satellite_scratch/output_data/contigs
cp ${files}/${sample_name}/${sample_name}.contigs.min500.fasta.gz /scratch/syersin2/Satellite_scratch/output_data/contigs
cp ${files}/${sample_name}/${sample_name}.contigs.min1000.fasta.gz /scratch/syersin2/Satellite_scratch/output_data/contigs

mkdir -p /scratch/syersin2/Satellite_scratch/output_data/scaffolds
cp ${files}/${sample_name}/${sample_name}.scaffolds.min0.fasta.gz /scratch/syersin2/Satellite_scratch/output_data/scaffolds
cp ${files}/${sample_name}/${sample_name}.scaffolds.min500.fasta.gz /scratch/syersin2/Satellite_scratch/output_data/scaffolds
cp ${files}/${sample_name}/${sample_name}.scaffolds.min1000.fasta.gz /scratch/syersin2/Satellite_scratch/output_data/scaffolds