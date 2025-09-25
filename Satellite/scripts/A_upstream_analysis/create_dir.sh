#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name create_dir
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G
#SBATCH --time 00:30:00

cd /scratch/syersin2/Satellite_scratch/data

for prefix in $(ls *.fastq.gz | sed -E 's/_FASTP_BOWTIE2_R[12]_001[.]fastq.gz//' | uniq)
	do
    mkdir ${prefix}
	mkdir -p /users/syersin2/Satellite/output/A_upstream_output/${prefix}
    cp ${prefix}_FASTP_BOWTIE2_R1_001.fastq.gz ${prefix}
    cp ${prefix}_FASTP_BOWTIE2_R2_001.fastq.gz ${prefix}
done