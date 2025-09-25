#!/bin/bash

#SBATCH --partition cpu
#SBATCH --job-name copy_files
#SBATCH --output /scratch/syersin2/Satellite_scratch/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G
#SBATCH --time 00:30:00

# Step done in terminal
# cp -r /nas/FAC/FBM/DMF/pvonaesc/default/D2c/Simon\ Yersin/Satellite/PacBio_Strept_Isolates/Fly_assemblies /scratch/syersin2/Satellite_scratch/Isolates/
# conda activate /work/FAC/FBM/DMF/pvonaesc/vonasch_lab_general/syersin/FastANI/fastani
#ls > genomes_list.txt
# fastANI --ql genomes_list.txt --rl genomes_list.txt -t 32 -o fastANI_PacBio.txt

## Variables
assemblies=/scratch/syersin2/Satellite_scratch/Isolates/Fly_assemblies
outdir=/scratch/syersin2/Satellite_scratch/Isolates/assemblies

cd ${assemblies}

for files in $(ls)
    do
    cd ${assemblies}/${files}
    cp assembly.fasta ${outdir}
    cd ${outdir}
    rename assembly ${files}_assembly assembly.fasta
    cd ${assemblies}
done

config=/users/syersin2/Satellite/scripts/Isolates/name_config.txt

cd ${outdir}
rename -v bc2158.seqwell_UDI3_Flye '' *.fasta

for files in $(ls *.fasta | sed -E 's/_assembly.fasta//')
    do
    sample_name=$(awk -v file_name="$files" '$4==file_name {print $3}' $config)
    rename -v "${files}" "${sample_name}" ${files}_assembly.fasta
done
