#!/bin/bash

#SBATCH --job-name ISP6
#SBATCH --output /scratch/syersin2/Satellite_scratch/transmission_dir/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/transmission_dir/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 150G
#SBATCH --time 10:00:00
#SBATCH --array=14,29,32,91,102,116

#Load modules
module load gcc/11.4.0
module load miniconda3/22.11.1

#Activate the Instrain environment
eval "$(/dcsrsoft/spack/external/micromamba/1.4.9/bin/micromamba-linux-64 shell hook --shell bash --root-prefix /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/mcreze/MAMBA/mamba_root 2> /dev/null)"
micromamba activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/mcreze/INSTRAIN/instrain

#InStrain profile takes as input a fasta file and a bam file (sorted and indexed! .bam.bai need to be next to the .bam) and runs a series of steps to characterize the nucleotide diversity, SNSs and SNVs, linkage, etc..
# If providing the .bam and .bam.bai, no need of samtools in inStrain 

# Define directories
indir=/scratch/syersin2/Satellite_scratch/transmission_dir/mapping_output
refdir=/scratch/syersin2/Satellite_scratch/transmission_dir/tables/genomes_db
filedir=/scratch/syersin2/Satellite_scratch/transmission_dir/tables
outdir=/scratch/syersin2/Satellite_scratch/transmission_dir/instrain/instrain_profile

# Change to input directory
cd ${indir}

# Define the array variables
sample=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)
rm -r ${outdir}/${sample}_inStrainOut

echo “###############################################”
echo $id “Start time: $(date -u) for sample: ” ${sample}
echo “##############################################”

inStrain profile ${indir}/${sample}/${sample}_sorted_allGenomesMapped.bam ${refdir}/allGenomes.fasta \
    -o ${outdir}/${sample}_inStrainOut \
    -s ${filedir}/genomes.stb \
    -g ${refdir}/prodigal/allGenomes.fna \
    -p 8 \
    --database_mode \
    --skip_plot_generation

echo “##############################################”
echo $id “Finishing time: $(date -u) ”
echo “##############################################”
