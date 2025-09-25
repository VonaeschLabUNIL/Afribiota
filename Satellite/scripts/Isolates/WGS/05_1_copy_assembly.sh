#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name copy_assemblies
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --time 00:05:00
#SBATCH --array=1-75

# Variables
indir=/scratch/syersin2/Satellite_scratch/Isolates/spades
outdir=/scratch/syersin2/Satellite_scratch/Isolates/wgs_assemblies
files_dir=/scratch/syersin2/Satellite_scratch/Isolates/spades_log

## Array variables
cd ${indir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Copy assemblies
cp ${indir}/${sample_name}/${sample_name}.scaffolds.min1000.fasta ${outdir}

## Copy files of interest 
mkdir -p ${files_dir}
cp ${indir}/${sample_name}/spades.log /${files_dir}/${sample_name}_spades.log
cp ${indir}/${sample_name}/scaffolds.fasta /${files_dir}/${sample_name}_scaffolds.fasta
cp ${indir}/${sample_name}/contigs.fasta /${files_dir}/${sample_name}_contigs.fasta
cp ${indir}/${sample_name}/${sample_name}.min200.assembly.stats /${files_dir}/
cp ${indir}/${sample_name}/${sample_name}.min1000.assembly.stats /${files_dir}/

## Zipping fasta files
gzip ${files_dir}/*.fasta