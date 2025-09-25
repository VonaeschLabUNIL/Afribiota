#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name backmapping_GADU
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --mail-type BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user simon.yersin@unil.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 200G
#SBATCH --time 12:00:00
#SBATCH --array=1-46

## Load module
module load gcc/11.4.0
module load bowtie2/2.5.1
module load samtools/1.17

## Variables
index_indir=/scratch/syersin2/Satellite_scratch/output_data/metaspades
reads_indir=/scratch/syersin2/Satellite_scratch/data
tmp=/scratch/syersin2/Satellite_scratch/tmp

## Array variables
config=/users/syersin2/Satellite/scripts/A_upstream_analysis/config_gastric_duodenum.txt
cd ${index_indir}
sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

## Backmapping
rm -r ${index_indir}/${sample_name}/${sample_name}_backmapping
mkdir -p ${index_indir}/${sample_name}/${sample_name}_backmapping

echo "#####################################################"
echo $id " Start time of the loop : $(date -u)"
echo "#####################################################"

## WARNING: +1 as awk reads the colname as row 1
for reads in $(awk 'NR >= 2 && NR <= 47 {print $2}' $config)
    do
    echo "Alignment of " ${reads} " on " ${sample_name} " is running and starts at $(date -u)" 
    bowtie2 -x ${index_indir}/${sample_name}/${sample_name}_index/${sample_name} \
        -1 ${reads_indir}/${reads}/${reads}_FASTP_BOWTIE2_R1_001.fastq.gz \
        -2 ${reads_indir}/${reads}/${reads}_FASTP_BOWTIE2_R2_001.fastq.gz \
        -S ${index_indir}/${sample_name}/${sample_name}_backmapping/${reads}_mapped_to_${sample_name}.sam \
        --threads 24
    echo "Changing alignment of " ${reads} " on " ${sample_name} " to sorted bam file at $(date -u)"
    samtools view -F 4 -bh ${index_indir}/${sample_name}/${sample_name}_backmapping/${reads}_mapped_to_${sample_name}.sam \
        | samtools sort -O bam -@ 4 -m 4G -T ${tmp} > ${index_indir}/${sample_name}/${sample_name}_backmapping/${reads}_mapped_to_${sample_name}.bam
    rm ${index_indir}/${sample_name}/${sample_name}_backmapping/${reads}_mapped_to_${sample_name}.sam
    echo "Alignment of " ${reads} " on " ${sample_name} " is completed at $(date -u)"
done

echo "#####################################################"
echo $id " Finishing time of the loop : $(date -u)"
echo "#####################################################"