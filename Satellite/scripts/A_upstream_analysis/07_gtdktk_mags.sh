#!/bin/bash

#SBATCH --job-name=GTDBTK
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time 4:00:00
#SBATCH --array=1-129

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment - DO NOT MODIFY!
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/GTDB_TK/gtdbtk

#### Create directories and define variables
gtdb_indir=/scratch/syersin2/Satellite_scratch/mags_output
keep_dir=/users/syersin2/Satellite

## Array variables
cd ${gtdb_indir}/MAGs
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

## Create output directory
mkdir -p ${gtdb_indir}/GTDB/${sample_name}

#run GTBD-Tk
gtdbtk classify_wf --genome_dir ${gtdb_indir}/MAGs/${sample_name} \
    --out_dir ${gtdb_indir}/GTDB/${sample_name} \
    --mash_db ${gtdb_indir}/GTDB/${sample_name}/gtdbtk.msh \
    --extension fa \
    --prefix ${sample_name}_gtdbtk \
    --cpus 32

mkdir -p /users/syersin2/Satellite/output/mags_analysis/GTDB
cp ${gtdb_indir}/GTDB/${sample_name}/${sample_name}_gtdbtk.ar53.summary.tsv ${keep_dir}/output/mags_analysis/GTDB
cp ${gtdb_indir}/GTDB/${sample_name}/${sample_name}_gtdbtk.bac120.summary.tsv ${keep_dir}/output/mags_analysis/GTDB
cp ${gtdb_indir}/GTDB/${sample_name}/${sample_name}_gtdbtk.log ${keep_dir}/Satellite/logs
cp ${gtdb_indir}/GTDB/${sample_name}/${sample_name}_gtdbtk.warnings.log ${keep_dir}/Satellite/logs