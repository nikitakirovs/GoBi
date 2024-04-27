import argparse
from ExtractBlastResults import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     'Extract potential CDS from a genome with tblastx results (has to be format 6)')
    parser.add_argument('-g', '--genome', required=True, help='Path to genome fasta file')
    parser.add_argument('-b', '--blast_results', required=True, help='Path to tblastx results in format 6')
    parser.add_argument('-o', '--output', required=True, help='File name for fasta file containing extracted CDS')
    args = parser.parse_args()

    blast_results = args.blast_results
    genome_path = args.genome

    tblastx_content = read_blast_results_to_list(blast_results)
    # get pyfaidx genome
    genome = Fasta(genome_path)

    records = get_tblastx_match_SeqIO_records(genome,tblastx_content)
    write_to_fasta(records, args.output)
