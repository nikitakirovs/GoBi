#!/bin/bash

# Check if the file argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <species_file>"
    exit 1
fi

# Read each species name from the file
while IFS= read -r species_name; do
    echo "Processing species: $species_name"

    ./commands2.sh "${species_name}"
done < "$1"

# write ./process_species_file.sh species_names.txt in order to download all the files for the species
