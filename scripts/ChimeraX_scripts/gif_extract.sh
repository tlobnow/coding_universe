#!/usr/bin/env bash

# The loop iterates over every GIF file in the directory.
# basename -- "$gif" .gif extracts the filename without the .gif extension.
# The ffmpeg command uses the -vf "select=eq(n\,0)" option to select the first frame (n=0) of the GIF.
# The resulting frame is saved as a PNG with the same name as the original GIF.

FILE_PATH=/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/PDB_BEST_iSCORE/BOTB

for gif in ${FILE_PATH}/*.gif; do
    # Extract the filename without the extension
    filename=$(basename -- "$gif" .gif)
    # Extract the first frame and save as PNG
    ffmpeg -i "$gif" -vf "select=eq(n\,0)" -q:v 3 "${FILE_PATH}/${filename}.png"
done
