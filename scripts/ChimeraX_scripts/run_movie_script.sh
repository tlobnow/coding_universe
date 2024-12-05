#!/usr/bin/env bash

#input_path=/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/PDB_BEST_iSCORE/BDLD_19_x6_model_5.pdb
input_path=/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/PDB_BEST_iSCORE/BOTB/BDLD_1_WP_013032525.1__Nitrosococcus_halophilus_x10_model_4.pdb
output_path="${input_path%.pdb}" # strips the .pdb suffix
temp_output="${output_path}.mp4"  # Temporary output for ChimeraX

# Execute ChimeraX with the script and pass the arguments
#chimeraX --script "./movie_script.cxc $input_path $temp_output"
chimeraX --script "./movie_script_x_only.cxc $input_path $temp_output"

# Convert the saved movie to GIF using ffmpeg
#ffmpeg -y -i "$temp_output" "${output_path}.gif"
ffmpeg -y -i "$temp_output" -vf "scale=320:-1" "${output_path}.gif"

# Optionally, remove the temporary mp4 file
rm "$temp_output"
