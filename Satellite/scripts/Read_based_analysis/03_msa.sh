#!/bin/bash

#SBATCH --job-name msa
#SBATCH --output /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/strainphlan/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 20G
#SBATCH --time 00:30:00
#SBATCH --array=1,3-481

# Module
module load gcc/11.4.0
module load miniconda3/22.11.1
module load samtools/1.17
module load blast-plus/2.14.1
module load mafft/7.505
module load raxml/8.2.12

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/MetaPhlan/metaphlan

## Variables
indir=/scratch/syersin2/Satellite_scratch/strainphlan/consensus_markers
outdir_clade=/scratch/syersin2/Satellite_scratch/strainphlan/CladeMarkers
outdir_sp=/scratch/syersin2/Satellite_scratch/strainphlan/StrainPhlAn_out

## Array variables
cd ${indir}
clade=$(awk -F'\t' 'NR>1 {print $1}' print_clades_only.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "Extracting markers: " ${clade}

# ## Extract the markers from MetaPhlAn database
extract_markers.py -c ${clade} \
     -o ${outdir_clade}/

# echo "Done"
echo "Build the multiple sequence alignment and the phylogenetic tree"

# ## Build the multiple sequence alignment and the phylogenetic tree
strainphlan -s ${indir}/*.json.bz2 \
    -m ${outdir_clade}/${clade}.fna \
    -o ${outdir_sp} \
    -c ${clade} \
    --phylophlan_mode accurate \
    --nproc 10

echo "Done for clade: " ${clade}