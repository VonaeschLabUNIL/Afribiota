#!/bin/bash

#SBATCH --job-name rename_contigs
#SBATCH --output /scratch/syersin2/Satellite_scratch/transmission_dir/std_output/%x_%j.out
#SBATCH --error /scratch/syersin2/Satellite_scratch/transmission_dir/std_output/%x_%j.err
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 00:10:00

cd /scratch/syersin2/Satellite_scratch/transmission_dir/dereplication/dereplicated_genomes

#!/bin/bash

# Loop through .fa files meeting the criteria
for file in *.fa; do
    # Process files starting with AFB or Ssal, or ending with assembly.fa
    if [[ "$file" == AFB* || "$file" == Ssal* || "$file" == *assembly.fa ]]; then
        # Extract the isolate name from the file name (up to the first period or underscore)
        isolate_name=$(basename "$file" | sed 's/\..*//; s/_.*//')

        # Create a temporary file to safely overwrite the original
        tmp_file=$(mktemp)

        # Add isolate name as a prefix to the headers
        awk -v isolate_name="${isolate_name}_" '/^>/ {print ">" isolate_name substr($0, 2); next} {print}' "$file" > "$tmp_file"

        # Replace the original file with the updated one
        mv "$tmp_file" "$file"

        echo "Processed and updated $file"
    fi
done
