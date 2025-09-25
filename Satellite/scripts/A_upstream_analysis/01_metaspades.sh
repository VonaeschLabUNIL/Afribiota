#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name metaSPAdes
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 400G
#SBATCH --time 48:00:00
#SBATCH --array=1-129

## Load modules
module load gcc/11.4.0
module load spades/3.15.5
module load python/3.10.13

## Variables
spades_workdir=/scratch/syersin2/Satellite_scratch
spades_indir=/scratch/syersin2/Satellite_scratch/data
spades_outdir=/scratch/syersin2/Satellite_scratch/output_data/metaspades
keep_dir=/users/syersin2/Satellite/output/A_upstream_output

## Array variables
cd ${spades_indir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)
cd ${spades_workdir}

## Create directory
rm -r ${spades_workdir}/output_data/metaspades/${sample_name}
mkdir -p ${spades_outdir}/${sample_name}

## Run metaSPAdes
spades.py --meta \
    --pe1-1 ${spades_indir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
    --pe1-2 ${spades_indir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R2_001.fastq.gz \
    -o ${spades_outdir}/${sample_name} \
    -k 21,33,55,77,99,127 \
    -t 32 \
    -m 400

# Cleaning files
rm -r ${spades_outdir}/${sample_name}/misc
rm -r ${spades_outdir}/${sample_name}/K*
rm -r ${spades_outdir}/${sample_name}/tmp

#gzip ${keep_dir}/${sample_name}/*.fasta
#gzip ${keep_dir}/${sample_name}/*.fastg