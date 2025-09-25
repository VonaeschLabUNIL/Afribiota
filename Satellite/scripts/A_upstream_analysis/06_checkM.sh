#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name CheckM_MAGs
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 36
#SBATCH --mem 250G
#SBATCH --time 48:00:00
#SBATCH --array=1-129

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/checkM/checkm

# Variables
qc_indir=/scratch/syersin2/Satellite_scratch/mags_output

## Array variables
cd ${qc_indir}/MAGs
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Run checkm
mkdir -p ${qc_indir}/checkM/${sample_name}

checkm lineage_wf ${qc_indir}/MAGs/${sample_name} \
    ${qc_indir}/checkM/${sample_name} \
    -x fa \
    -t 36 \
    -f ${qc_indir}/checkM/${sample_name}/${sample_name}_checkM_stats.tsv \
    --tab_table

# Saving files 
cd ${qc_indir}/checkM/${sample_name}
rename -v checkm "${sample_name}" checkm.log
mkdir -p /users/syersin2/Satellite/output/mags_analysis/checkM
cp ${qc_indir}/checkM/${sample_name}/${sample_name}_checkM_stats.tsv /users/syersin2/Satellite/output/mags_analysis/checkM
cp ${qc_indir}/checkM/${sample_name}/${sample_name}.log /users/syersin2/Satellite/logs