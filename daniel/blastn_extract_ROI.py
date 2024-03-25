import argparse
from ExtractBlastResults import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     'Extract Regions Of Interest from a genome with blastn results (has to be format 6)')
    parser.add_argument('-g', '--genome', required=True, help='Path to genome fasta file')
    parser.add_argument('-b', '--blast_results', required=True, help='Path to blastn results in format 6')
    parser.add_argument('-o', '--output', required=True, help='File name for fasta file containing extracted ROIs')
    args = parser.parse_args()

    blast_results = args.blast_results
    genome_path = args.genome

    blastn_content = read_blast_results_to_list(blast_results)
    unique_landmarks = get_unique_landmarks(blastn_content)
    genome = Fasta(genome_path)
    records = get_landmark_SeqIO_records(genome, unique_landmarks)
    write_to_fasta(records, args.output)
