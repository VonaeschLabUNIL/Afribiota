#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name pairwisedist
#SBATCH --output /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 00:02:00
#SBATCH --array 1-480 # How many SGB do we have?

# Module
#module load gcc/11.4.0
module load python/3.8.18

# Activate venv environment
source /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/PyPhlAn/venv/bin/activate

## Variables
indir=/scratch/syersin2/Satellite_scratch/strainphlan/StrainPhlAn_out
outdir=/scratch/syersin2/Satellite_scratch/strainphlan/Phylo_dist
pythondir=/work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/PyPhlAn

mkdir -p ${outdir}

## Define array variable
# Extract the filename corresponding to this SLURM array task ID
cd ${indir}
file=$(ls RAxML_bestTree.t__SGB*.StrainPhlAn4.tre | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract the unique identifier from the filename
SGB=$(basename "${file}" | sed 's/^RAxML_bestTree\.//' | sed 's/\.StrainPhlAn4\.tre$//')

# Debug: Print the identifier
echo “Processing identifier: $SGB”

## Pairwise phylogenetic distances
python ${pythondir}/tree_pairwisedists.py -n RAxML_bestTree.${SGB}.StrainPhlAn4.tre ${outdir}/${SGB}_nGD.tsv

echo "Done identifying:$SGB"