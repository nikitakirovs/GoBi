from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def parse_gff(gff_file):
    cds_features = []
    with open(gff_file) as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'CDS':
                    start, end = int(parts[3]), int(parts[4])
                    strand = parts[6]
                    attributes = parts[8].split(';')
                    # Extracting the ID from the attributes field
                    cds_id = [x.split('=')[1] for x in attributes if x.startswith('ID=')][0]
                    cds_features.append((cds_id, (start, end, strand)))
    return cds_features

# concat genome to one string
def load_genome(genome_path):
    with open(genome_path) as f:
        genome_record = SeqIO.parse(genome_path, 'fasta')
    return genome_record

# take cds_features list from parse_gff and return a list only containing the cds from the specified range
def cds_from_range(range_start, range_end, cds_features):
    cds_to_return = []
    for tpl in cds_features:
        id = tpl[0]
        cds_start, cds_end, _ = tpl[1]
        if (range_start <= cds_start <= cds_end and range_start <= cds_end <= range_end ):
            cds_to_return.append(tpl)

    return cds_to_return

def extract_cds_sequences(genome, cds_list):
    cds_sequences = []
    for i, location in cds_list:
        start, end, strand = location
        if strand == '+':
            cds_seq = genome[start - 1:end].seq
        else:
            cds_seq = genome[start - 1:end].seq.reverse_complement()
        cds_record = SeqRecord(cds_seq, id=f'{i}', description=f'CDS_{start}_{end}')
        cds_sequences.append(cds_record)
    return cds_sequences


if __name__ == '__main__':
    # specify path of gff and genome
    gff_path = '/home/daniel/blast/test-data/apis-mellifera/genomic.gff'
    genome_path = '/home/daniel/blast/test-data/apis-mellifera/GCF_003254395.2_Amel_HAv3.1_genomic.fna'

    # set range of interest. a range since we look at the sequence of yellow-e3 to MRJP9
    yellowe3_start = 2259179
    mrjp9_end = 2325880
    flank_size = 50000

    # extract all cds from gff
    cds_features = parse_gff(gff_path)
    print(f'CDS count whole genome: {len(cds_features)}')

    # reduce to cds in range of interest
    cds_filtered = cds_from_range(yellowe3_start-flank_size,mrjp9_end+flank_size,cds_features)
    print(f'CDS in range of interest: {len(cds_filtered)}')

    # load genome
    genome = load_genome(genome_path)

    # extract cds from genome
    cds_records = extract_cds_sequences(genome,cds_filtered)

    #write cds to file
    output_file = 'apis-mellifera-cds-selection.fasta'
    SeqIO.write(cds_records, output_file, 'fasta')

    print(f'CDS of interest written to file: {output_file}')
