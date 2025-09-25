#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --job-name=GTDBTK_ISO
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --time 4:00:00
#SBATCH --export=None

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1

# Activate conda environment - DO NOT MODIFY!
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/GTDB_TK/gtdbtk

#### Create directories and define variables
inputfolder=/scratch/syersin2/Satellite_scratch/Isolates/assemblies

#Create output directory and make sure is empty
outdir=/scratch/syersin2/Satellite_scratch/Isolates/outputs/GTDB
rm -r ${outdir} 
mkdir ${outdir}

#run GTBD-Tk
gtdbtk classify_wf --genome_dir ${inputfolder} --out_dir ${outdir} --mash_db ${outdir}/gtdbtk.msh --extension fasta --cpus 16