#!/usr/bin/env bash

echo "Starting Genome Downloads"

GENOME_DIR="/Users/u_lobnow/Documents/Github/coding_universe/SOTA/01_raw_data/GENOMES/"
DOWNLOAD_LIST="${GENOME_DIR}/download_list.csv"

echo "Downloading Data from $DOWNLOAD_LIST"

while IFS=$'\t' read -r assembly; do
	
	ZIP_FILE="${GENOME_DIR}/${assembly}.zip"
	EXTRACTED_FOLDER="${GENOME_DIR}/${assembly}_extracted"

	if [ -f "$ZIP_FILE" ] || [ -d "$EXTRACTED_FOLDER" ]; then
        	echo "Assembly $assembly already downloaded and/or extracted, skipping..."
        	continue
    	fi

	echo "Downloading genome with assembly number $assembly"

	datasets download genome accession "$assembly" --include gff3,genome --filename "${GENOME_DIR}/${assembly}.zip"

	if [ -f "$ZIP_FILE" ]; then
		echo "Unzipping ${assembly}.zip"
		unzip -o "$ZIP_FILE" -d "$EXTRACTED_FOLDER"
	else
		echo "Unzipping failed for assembly $assembly"
	fi

done < $DOWNLOAD_LIST

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
