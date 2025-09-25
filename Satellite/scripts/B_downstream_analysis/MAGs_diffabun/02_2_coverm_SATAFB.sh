#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name coverm
#SBATCH --output /scratch/syersin2/SATAFB_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/SATAFB_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 50G
#SBATCH --time 06:00:00
#SBATCH --array=60-62

# Module
module load miniforge3/4.8.3-4-Linux-x86_64

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/CoverM/coverm

# Variables
# Strains
genome_dir=/scratch/syersin2/SATAFB_scratch/drep_strains/dereplicated_genomes
index_dir=/scratch/syersin2/SATAFB_scratch/drep_strains/index
reads_dir=/scratch/syersin2/SATAFB_scratch/reads
outdir=/scratch/syersin2/SATAFB_scratch/coverm_strains/coverm_out
TMPDIR=/scratch/syersin2/SATAFB_scratch/coverm_strains/tmp
# Species
#genome_dir=/scratch/syersin2/SATAFB_scratch/drep_species/dereplicated_genomes
#index_dir=/scratch/syersin2/SATAFB_scratch/drep_species/index
#reads_dir=/scratch/syersin2/Afribiota_scratch/reads
#outdir=/scratch/syersin2/SATAFB_scratch/coverm_species/coverm_out
#TMPDIR=/scratch/syersin2/SATAFB_scratch/coverm_species/tmp

## Array variables
cd ${reads_dir}
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p ${outdir}/${sample_name}

# Strains
coverm genome --coupled ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R1_001.fastq.gz \
     ${reads_dir}/${sample_name}/${sample_name}_FASTP_BOWTIE2_R2_001.fastq.gz \
     --genome-fasta-directory ${genome_dir} \
     -r ${index_dir}/SAT_AFB_all_genomes.fa \
     -x fa \
     -m tpm \
     -p bwa-mem \
     -t 16 \
     -o ${outdir}/${sample_name}/${sample_name}_CoverM_SATAFB.tsv

# Species
#coverm genome --coupled ${reads_dir}/${sample_name}/${sample_name}_bowtie2_R1_001.fastq.gz \
#     ${reads_dir}/${sample_name}/${sample_name}_bowtie2_R2_001.fastq.gz \
#     --genome-fasta-directory ${genome_dir} \
#     -r ${index_dir}/SAT_AFB_drep_species.fa \
#     -x fa \
#     -m tpm \
#     -p bwa-mem \
#     -t 16 \
#     -o ${outdir}/${sample_name}/${sample_name}_CoverM_SATAFB_species.tsv

