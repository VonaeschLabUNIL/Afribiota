#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name covermsp
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/mags_output/coverm_species/log/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 50G
#SBATCH --time 08:00:00
#SBATCH --array=1-129

# Module
# module load gcc/12.3.0
module load miniforge3/4.8.3-4-Linux-x86_64

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/CoverM/coverm

# Variables
genome_dir=/scratch/syersin2/Satellite_scratch/mags_output/dRep/dereplicated_genomes
index_dir=/scratch/syersin2/Satellite_scratch/mags_output/dRep/index
reads_dir=/scratch/syersin2/SATAFB_scratch/reads
outdir=/scratch/syersin2/Satellite_scratch/mags_output/coverm_species
TMPDIR=/scratch/syersin2/Satellite_scratch/mags_output/coverm_species/tmp

## Array variables
sample_name=$(ls ${reads_dir} | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p ${outdir}/${sample_name}

# Coverage:
# potential to use covered_fraction and add --min-covered-fraction 0
# Redo at species level
coverm genome --coupled ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
    ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R2_001.fastq.gz \
    --genome-fasta-directory ${genome_dir} \
    -r ${index_dir}/Satellite_drep_species.fa \
    -x fa \
    -m tpm \
    -p bwa-mem \
    -t 16 \
    -o ${outdir}/${sample_name}/${sample_name}_CoverM_tpm_species.tsv

    # Coverage:
# # potential to use covered_fraction and add --min-covered-fraction 0
# coverm genome --coupled ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
#     ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R2_001.fastq.gz \
#     --genome-fasta-directory ${genome_dir}/${sample_name} \
#     -x fa.gz \
#     -m covered_fraction \
#     --min-covered-fraction 0 \
#     -p bwa-mem \
#     -t 16 \
#     -o ${outdir}/${sample_name}/${sample_name}_CoverM_coverfract.tsv

# # Count might be necessary:
# # coverm genome --methods count --min-covered-fraction 0 --coupled
# coverm genome --coupled ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
#      ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R2_001.fastq.gz \
#      --genome-fasta-directory ${genome_dir}/${sample_name} \
#      -x fa.gz \
#      -m count \
#      --min-covered-fraction 0 \
#      -t 16 \
#      -o ${outdir}/${sample_name}/${sample_name}_CoverM_count.tsv
