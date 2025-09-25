#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name zipping
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 8G
#SBATCH --time 4:00:00
#SBATCH --array=1-129

## Variables
files=/scratch/syersin2/Satellite_scratch/output_data/metaspades

## Array variables
cd ${files}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

#cd ${files}/${sample_name}/
#gzip *.fasta
#gzip *.gfa
#gzip *.fastg

mags_file=/scratch/syersin2/Satellite_scratch/mags_output/MAGs
cd ${mags_file}/${sample_name}/
gzip *.fa