# Repository of the paper "Exploring the phylogeny of the mrjps and the yellow gene family"

## Extraction of species

Run the script process_species_file.sh as following: 

```shell

./process_species_file.sh species_names.txt

```
The script expects a species_list as input.  
The script creates directories with data files, that where found by searching with the species name and taking the first entry with the NCBI Datasets command-line tool. Therefore the package has to be installed in the environment. (see https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)  
  
#### Restriction of species list: 

The species_list is a .txt file containing the complete Latin names of all species analyzed in the paper. Each name is separated by a newline (\n).  

Example: 

```shell
Apis mellifera
Drosophila melanogaster
Bombyx mori
etc. 

```



## Running the blast pipeline

```shell

python pipeline.py -g <path_subdirectory> -s <path_to_species.txt> -d <output_dir_path> -a <dir_alias_db> -m <blosum_mat>

```
Explanation to options: 
- -g: expects path to directory containing one subdirectory for each species with its genomic.fna  
- -s: expects path to the species.txt  
- -d: directory path to save database for each species, used to create an alias database linking them together  
- -a: directory path for alias database  
- -m: int n specifying BLOSUM n matrix for tblastx  
