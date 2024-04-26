#!/bin/bash

# Check if the parameter is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <species_name>"
    exit 1
fi

species_name="$1"

# Run the command with the provided species name and store the output
summary_output=$(datasets summary genome taxon "$species_name")

# Extract the accession ID from the summary output
accession_id=$(echo "$summary_output" | grep -o '"accession":"[^"]*' | head -n 1 | cut -d ':' -f 2| cut -c 2-)

if [ -n "$accession_id" ]; then
    echo "Accession ID found: $accession_id"

    # Download using the extracted accession ID
    output_file="${1}_cds.fna"
    datasets download genome accession "$accession_id" --include cds
    unzip -j ncbi_dataset.zip "*/cds_from_genomic.fna" -d "$(dirname "$0")"
    mv "$(dirname "$0")/cds_from_genomic.fna" "$(dirname "$0")/$output_file"
    rm ncbi_dataset.zip
    echo "Genome data downloaded and saved as: $output_file"
else
    echo "Accession ID not found"
fi

# in order to run the script write: 
# ./script_name.sh "Species name"


