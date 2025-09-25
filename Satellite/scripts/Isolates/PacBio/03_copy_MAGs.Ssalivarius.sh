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
assemblies=/scratch/syersin2/Satellite_scratch/transmission_dir/MAGs
outdir=/scratch/syersin2/Satellite_scratch/Isolates/assemblies
config=/users/syersin2/Satellite/scripts/Isolates/S.salivarius_MAGs_dRep.txt

cd ${assemblies}

for files in $(awk -v file_name="$files" '{print $1}' $config)
    do
    cp ${files} ${outdir}
done
