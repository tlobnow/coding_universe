#!/usr/bin/env bash

echo "Starting Rename Protocol"

GENOME_DIR="/Users/u_lobnow/Documents/Github/coding_universe/SOTA/01_raw_data/GENOMES/"
DOWNLOAD_LIST="${GENOME_DIR}/download_list.csv"

# Loop through each *_extracted directory
for extracted_dir in "$GENOME_DIR"/*_extracted; do
    if [ -d "$extracted_dir" ]; then
        echo "Processing $extracted_dir"

        # Extract the assembly prefix from the folder name
        assembly_prefix=$(basename "$extracted_dir" | sed 's/_extracted$//')
	echo $assembly_prefix

        # Find all files in the extracted_dir, excluding the .fna file
        #find "$extracted_dir" -maxdepth 1 -type f ! -name "\*.fna" | while read -r file; do
	find "$extracted_dir" -maxdepth 1 -type f ! -name "${assembly_prefix}_*" | while read -r file; do
            filename=$(basename "$file")
            new_filename="${assembly_prefix}_${filename}"
            mv "$file" "${extracted_dir}/${new_filename}"
            echo "Renamed $filename to $new_filename"
        done
    fi
done

echo "All files renamed"
