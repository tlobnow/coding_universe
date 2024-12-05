#!/usr/bin/env bash

# Specify the directory containing the PDB files
#directory_path="/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/PDB_BEST_iSCORE/BOTB/"
#directory_path="/Users/u_lobnow/Documents/Github/transferGit/MYD88_MOUSE_x1_P68373_x1/MODELS/"
directory_path="/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/MYD88/MYD88_MOUSE_x1_A2A3V1_x1_rep1"

# Loop over all PDB files in the directory
for input_path in "${directory_path}/"*.pdb; do
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
done
