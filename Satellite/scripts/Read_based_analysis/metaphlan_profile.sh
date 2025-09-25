#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name metaphlan
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 20G
#SBATCH --time 4:00:00
#SBATCH --array=1-129

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
mpa_workdir=/scratch/syersin2/Satellite_scratch
mpa_indir=/scratch/syersin2/Satellite_scratch/data
## Array variables
cd ${mpa_indir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir ${mpa_workdir}/output_data/METAPHLAN

echo "Performing metaphlan3 analysis on sample": ${sample_name}
metaphlan ${mpa_indir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
    --nproc 16 \
    --no_map \
    -t rel_ab_w_read_stats \
    --input_type fastq \
    -o ${mpa_workdir}/output_data/METAPHLAN/${sample_name}_metaphlan_profile.txt

## Change to GTDB taxonomy
mkdir ${mpa_workdir}/output_data/METAPHLAN/GTDB_tax
sgb_to_gtdb_profile.py -i ${mpa_workdir}/output_data/METAPHLAN/${sample_name}_metaphlan_profile.txt \
    -o ${mpa_workdir}/output_data/METAPHLAN/GTDB_tax/${sample_name}_gtdb.txt
