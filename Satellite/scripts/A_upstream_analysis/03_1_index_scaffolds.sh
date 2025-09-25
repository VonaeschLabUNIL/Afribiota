#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name build_index
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 20G
#SBATCH --time 04:00:00
#SBATCH --array=1-129

## Load module
module load gcc/11.4.0
module load bowtie2/2.5.1

## Variables
index_indir=/scratch/syersin2/Satellite_scratch/output_data/metaspades

## Array variables
cd ${index_indir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

## Build index
mkdir -p ${index_indir}/${sample_name}/${sample_name}_index
bowtie2-build --threads 16 ${index_indir}/${sample_name}/${sample_name}.scaffolds.min1000.fasta ${index_indir}/${sample_name}/${sample_name}_index/${sample_name}

#mkdir -p /users/syersin2/Satellite/output/A_upstream_output/${sample_name}/${sample_name}_index
#cp ${index_indir}/${sample_name}/${sample_name}_index/*.bt2 /users/syersin2/Satellite/output/A_upstream_output/${sample_name}/${sample_name}_index/