#!/usr/bin/env bash

GENOME_DIR="/Users/u_lobnow/Documents/Github/coding_universe/SOTA/01_raw_data/GENOMES/"
DOWNLOAD_LIST="${GENOME_DIR}/download_list.csv"

echo "Flattening folder structure in each extracted genome folder"

# Loop through each *_extracted directory
for extracted_dir in "$GENOME_DIR"/*_extracted; do
    if [ -d "$extracted_dir" ]; then

        # Find all files in the extracted_dir, excluding directories, and move them to the extracted_dir
	find "$extracted_dir" -mindepth 2 -type f | while read -r file; do mv "$file" "$extracted_dir" ; done

        # Remove empty directories
        find "$extracted_dir" -type d -empty -delete
    fi
done

echo "All folder structures flattened"
