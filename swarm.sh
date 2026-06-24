#!/bin/bash

#######################
### Setup
###

INPUT_DIR="10000_species_10_dynamic/rsearch_results/fasta"  # Path to the directory containing input FASTA files
FILE_EXT="*.fasta"                                          # File extension pattern used to match input FASTA files
OUTPUT_DIR="10000_species_10_dynamic/swarm_results"         # Path to the directory where all results will be stored

THREADS=4                                                   # Number of CPU cores to use for parallel processing

SWARM_D=1                                                   # Maximum number of differences allowed between sequences within a swarm (default: 1)

mkdir -p "$OUTPUT_DIR"                                      # Create the output directory if it does not already exist

#######################
### Main Loop
###

for INPUT_FASTA in "$INPUT_DIR"/$FILE_EXT; do           # Iterate over all FASTA files in the input directory
  
  BASENAME=$(basename "$INPUT_FASTA" .fasta)            # Extract the filename without path or extension, used as a base for output filenames
  UC_FILE="${OUTPUT_DIR}/${BASENAME}_swarm.uc"          # Define the path for the SWARM UC (uclust-format) output file

  echo "=========================================="
  echo " Starter pipeline for: $INPUT_FASTA"
  echo " Resultater lagres i mappen: $OUTPUT_DIR/"
  echo " Tråder: $THREADS | Swarm d-verdi: $SWARM_D"
  echo "=========================================="
  
  # Clustering with Swarm
  swarm -d "$SWARM_D" \                                    # Set the maximum number of differences allowed between sequences in a swarm
        -f \                                               # Enable the fastidious option to reduce the number of small spurious OTUs
        -z \                                               # Use USEARCH/VSEARCH-compatible output format for the UC file
        -t "$THREADS" \                                    # Set the number of threads for parallel execution
        -u "$UC_FILE" \                                    # Write cluster assignments to the UC file
        "$INPUT_FASTA"                                     # Input FASTA file to cluster

  # Error handling
  if [ $? -ne 0 ]; then                                                                 # Check the exit code of the previous command; non-zero indicates an error
      echo "   [ERROR] Something went wrong with $INPUT_FASTA. Skipping to next file."
      continue                                                                          # Skip to the next file in the loop instead of terminating the entire script
  fi

  echo "   [OK] Saved as $BASENAME"
  
done