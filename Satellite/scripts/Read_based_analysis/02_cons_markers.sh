#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name SP_CM
#SBATCH --output /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G
#SBATCH --time 00:20:00
#SBATCH --array=91

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1
module load samtools/1.17

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
workdir=/scratch/syersin2/Satellite_scratch/strainphlan
indir=/scratch/syersin2/Satellite_scratch/strainphlan/sams
outdir=/scratch/syersin2/Satellite_scratch/strainphlan/consensus_markers

#rm -r ${outdir}
#mkdir -p ${outdir}

## Array variables
cd ${workdir}/Cleaned_reads
sample_name=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Produce consensus marker files for sample: " ${sample_name}

sample2markers.py -i ${indir}/${sample_name}.sam.bz2 \
     -o ${outdir} \
     -n 4

echo "Done for samples: " ${sample_name}

# Run then in terminal in the outdir:
# activate metaphlan environment
# strainphlan -s ./*.json.bz2 --print_clades_only -o . --nproc 8 > ./clades.txt
