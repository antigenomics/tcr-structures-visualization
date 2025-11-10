#!/bin/bash

# Usage: ./script.sh <source_directory> <output_directory>
# source_directory: The directory containing the folders to iterate through
# output_directory: The directory where the renamed .pdb files will be copied

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <source_directory> <output_directory>"
    exit 1
fi

SOURCE_DIR="$1"
OUTPUT_DIR="$2"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Iterate through each item in the source directory
for folder_path in "$SOURCE_DIR"/*; do
    if [ -d "$folder_path" ]; then
        folder=$(basename "$folder_path")
        subfolder="${folder_path}/${folder}_pmhc_oc"
        
        if [ -d "$subfolder" ]; then
            pdb_file="${subfolder}/ranked_0.pdb"
            
            if [ -f "$pdb_file" ]; then
                cp "$pdb_file" "${OUTPUT_DIR}/${folder}.pdb"
                echo "Copied ${pdb_file} to ${OUTPUT_DIR}/${folder}.pdb"
            else
                echo "File ranked_0.pdb not found in ${subfolder}"
            fi
        else
            echo "Subfolder ${folder}_pmhc_oc not found in ${folder_path}"
        fi
    fi
done