#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --job-name=orthofinder
#SBATCH --output /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/Isolates/std_output/%x_%j.err
#SBATCH --time 2:00:00
#SBATCH --export=None

# Module
module load gcc
module load miniconda3

# ### Create directories and define variables
bins=/scratch/syersin2/Satellite_scratch/Isolates/phylogeny/genomesFAA
workdir=/scratch/syersin2/Satellite_scratch/Isolates/phylogeny
rm -r ${workdir}/OrthoFinder
mkdir ${workdir}

echo "#####################################################"
echo $id "Starting time: $(date -u)"
echo "#####################################################"

################# Running Orthofinder #################
# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/OrthoFinder/orthofinder

orthofinder -f ${bins} -t 24 -a 24 -o ${workdir}/OrthoFinder -M msa

echo "#####################################################"
echo $id " Finishing time : $(date -u)"
echo "#####################################################"