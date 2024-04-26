# Repository of the paper "Exploring the phylogeny of the mrjps and the yellow gene family"

**Extraction of species**

Run the script process_species_file.sh as following: 

```shell

./process_species_file.sh species_names.txt

```

species_name.txt represents the list of species. #### here change! 



**Running the blast pipeline**

```shell
'python pipeline.py -g <Path to directory containing one subdirectory for each species with its genomic.fna> -s <path to species.txt> -d <dir to save db for each species to, used to create alias db linking them together> -a <dir for alias db> -m <int n specifying BLOSUM n matrix for tblastx>


```
