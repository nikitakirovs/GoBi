import os.path
import argparse
import subprocess

def file_is_empty(file_path):
    return os.path.getsize(file_path) == 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     'Do all the stuff')
    parser.add_argument('-g', '--genomes_dir', required=True,
                        help='Path to directory containing one subdirectory for each species with its genomic.fna')
    parser.add_argument('-s', '--species_file', required=True,
                        help='Unix ASCII file containing one species name on each line')
    parser.add_argument('-d', '--dbs_dir', required=True,
                        help='Blast db with subject CDS')
    parser.add_argument('-a', '--aliasdb_dir', required=True,
                        help='aliasdb dir')
    parser.add_argument('-m', '--matrix', required=True,
                        help='specifies int X for BLOSUM{X} matrix to be used')
    args = parser.parse_args()

    current_db = args.aliasdb_dir

    if os.path.exists(args.genomes_dir) and os.path.isdir(args.genomes_dir):
        genomes_dir = args.genomes_dir
    else:
        ValueError(f'{args.genomes_dir} is not a valid directory')

    if os.path.exists(args.species_file) and os.path.isfile(args.species_file):
        species_file = args.species_file
    else:
        ValueError(f'{args.species_file} is not an existing file')

    if os.path.exists(args.dbs_dir) and os.path.isdir(args.dbs_dir):
        dbs_dir = args.dbs_dir
    else:
        ValueError(f'{args.dbs_dir} does not exist or isn\'t a directory')

    matrix = f'BLOSUM{args.matrix}'

    species_list = []
    with open(species_file) as f:
        for line in f.readlines():
            species_list.append(line.strip())


    for species in species_list:
        species_dir = species.replace(' ', '-')
        current_genome_dir = f'{genomes_dir}/{species_dir}'
        genome_file = 'genome.fna'
        current_genome = f'{current_genome_dir}/{genome_file}'

        blastn_cmd = f'blastn -db {current_db}/aliasdb -query {current_genome} -gapopen 6 -gapextend 3 -penalty -3 -outfmt 11 -out {current_genome_dir}/{species_dir}_blastn_results-archive -num_threads 8'
        blastn_fmt6 = f'blast_formatter -archive {current_genome_dir}/{species_dir}_blastn_results-archive -outfmt 6 > {current_genome_dir}/{species_dir}_blastn_results'
        blastn_fmt5 = f'blast_formatter -archive {current_genome_dir}/{species_dir}_blastn_results-archive -outfmt 5 > {current_genome_dir}/{species_dir}_blastn_results.xml'
        extract_blastn = f'python blastn_extract_ROI.py -g {current_genome} -b {current_genome_dir}/{species_dir}_blastn_results -o {current_genome_dir}/{species_dir}_ROI.fasta'
        tblastx_cmd = f'tblastx -db {current_db}/aliasdb -query {current_genome_dir}/{species_dir}_ROI.fasta -outfmt 11 -matrix {matrix} -out {current_genome_dir}/{species_dir}_tblastx_results-archive -num_threads 8'
        tblastx_fmt6 = f'blast_formatter -archive {current_genome_dir}/{species_dir}_tblastx_results-archive -outfmt 6 > {current_genome_dir}/{species_dir}_tblastx_results'
        tblastx_fmt5 = f'blast_formatter -archive {current_genome_dir}/{species_dir}_tblastx_results-archive -outfmt 5 > {current_genome_dir}/{species_dir}_tblastx_results.xml'

        extract_tblastx = f'python tblastx_extract_seqs.py -g {current_genome} -b {current_genome_dir}/{species_dir}_tblastx_results -o {current_genome_dir}/{species_dir}_potential_CDS.fasta'
        make_dir = f'mkdir {dbs_dir}/{species_dir}'
        create_db = f'makeblastdb -in {current_genome_dir}/{species_dir}_potential_CDS.fasta -dbtype nucl -out {dbs_dir}/{species_dir}/{species_dir}'
        update_db = f'python update_current_db.py -d {dbs_dir} -a {current_db}'

        # blastn
        print(f'running: {blastn_cmd}')
        blastn_res = f'{current_genome_dir}/{species_dir}_blastn_results'
        print(subprocess.run(blastn_cmd, shell=True, capture_output=True).stderr.decode())

        print('Format archive output to format 6')
        print(f'running: {blastn_fmt6}')
        print(subprocess.run(blastn_fmt6, shell=True, capture_output=True).stderr.decode())

        print('Format archive output to format 5')
        print(f'running: {blastn_fmt5}')
        print(subprocess.run(blastn_fmt5, shell=True, capture_output=True).stderr.decode())
        # input('test')

            # check if blast results are empty
        if file_is_empty(blastn_res):
            print(f'blastn: No Significant alignments found in {species}')
            # input('Continue?')
            # print('Continuing with next species without updating aliasdb')
            continue



        print(f'running: {extract_blastn}')
        print(subprocess.run(extract_blastn, shell=True, capture_output=True).stderr.decode())

        # tblastx
        print(f'running: {tblastx_cmd}')
        tblastx_res = f'{current_genome_dir}/{species_dir}_tblastx_results' 
        print(subprocess.run(tblastx_cmd, shell=True, capture_output=True).stderr.decode())

        print('Format archive output to format 6')
        print(f'running: {tblastx_fmt6}')
        print(subprocess.run(tblastx_fmt6, shell=True, capture_output=True).stderr.decode())
        print('Format archive output to format 5')
        print(f'running: {tblastx_fmt5}')
        print(subprocess.run(tblastx_fmt5, shell=True, capture_output=True).stderr.decode())

            # check if blast results are empty
        if file_is_empty(tblastx_res):
            print(f'tblastx: No Significant alignments found in {species}')
            # input('Continue?')
            # print('Continuing with next species without updating aliasdb')
            continue




        print(f'running: {extract_tblastx}')
        print(subprocess.run(extract_tblastx, shell=True, capture_output=True).stderr.decode())
        # input('Press enter to continue')
        print(f'running: {make_dir}')
        print(subprocess.run(make_dir, shell=True, capture_output=True).stderr.decode())
        print(f'running: {create_db}')
        print(subprocess.run(create_db, shell=True, capture_output=True).stderr.decode())
        print(f'running: {update_db}')
        print(subprocess.run(update_db, shell=True, capture_output=True).stderr.decode())


    print('done!')
