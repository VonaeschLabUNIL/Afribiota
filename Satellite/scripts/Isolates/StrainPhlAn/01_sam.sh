#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name SP_metaphlan
#SBATCH --output /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 50G
#SBATCH --time 01:00:00
#SBATCH --array=1-129

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
indir=/scratch/syersin2/Satellite_scratch/Isolates/cleaned_reads
outdir=/scratch/syersin2/Satellite_scratch/Isolates/strainphlan


## Array variables
cd ${indir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p | sed 's/\.fastq\.gz$//')

echo "Performing metaphlan4 analysis on sample": ${sample_name}


metaphlan ${indir}/${sample_name}/${sample_name}_FASTP_R1_001.fastq.gz
    --nproc 8 \
    --input_type fastq \
    -s ${outdir}/sams/${sample_name}.sam.bz2 \
    --bowtie2out ${outdir}/bowtie2/${sample_name}.bowtie2.bz2 \
    -o ${outdir}/profiles/${sample_name}_profiled.tsv

echo "Done for sample": ${sample_name}

# Need 104 files at the end