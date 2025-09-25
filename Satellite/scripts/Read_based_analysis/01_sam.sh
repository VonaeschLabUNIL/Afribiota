#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name SP_metaphlan
#SBATCH --output /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 30G
#SBATCH --time 02:00:00
#SBATCH --array=14,21,39,72,91,102,116

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
indir=/scratch/syersin2/Satellite_scratch/strainphlan/Cleaned_reads
outdir=/scratch/syersin2/Satellite_scratch/strainphlan

## Array variables
cd ${indir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Performing metaphlan4 analysis on sample": ${sample_name}

metaphlan ${indir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
    --nproc 8 \
    --input_type fastq \
    -s ${outdir}/sams/${sample_name}.sam.bz2 \
    --bowtie2out ${outdir}/bowtie2/${sample_name}.bowtie2.bz2 \
    -o ${outdir}/profiles/${sample_name}_profiled.tsv

echo "Done for sample": ${sample_name}