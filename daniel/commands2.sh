#!/bin/bash

# Check if the parameter is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <species_name>"
    exit 1
fi

species_name=$1
# Run the command with the provided species name and store the output
summary_output=$(datasets summary genome taxon "$species_name")
# Extract the accession ID from the summary output
accession_id=$(echo $summary_output | jq '.reports.[0]["accession"]')
accession_id="${accession_id//\"/}"
echo $accession_id

if [ -n "$accession_id" ]; then
    echo "Accession ID found: $accession_id"
    species_dir=${species_name//\ /-}
    # Download using the extracted accession ID
    echo create directory $species_dir
    mkdir $species_dir
    cd $species_dir
    datasets download genome accession $accession_id --include genome
    unzip -j ncbi_dataset.zip
    rm ncbi_dataset.zip *json *jsonl README.md
    echo "$species_name Genome downloaded and saved to: $species_dir"
else
    echo "Accession ID not found"
fi

# in order to run the script write: 
# ./script_name.sh "Species name"


