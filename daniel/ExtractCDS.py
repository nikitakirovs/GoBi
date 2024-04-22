from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta


def parse_gff(gff_file,landmark):
    cds_features = []
    with open(gff_file) as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'CDS' and parts[0] == landmark:
                    start, end = int(parts[3]), int(parts[4])
                    strand = parts[6]
                    attributes = parts[8].split(';')
                    # Extracting the ID from the attributes field
                    cds_id = [x.split('=')[1] for x in attributes if x.startswith('ID=')][0]
                    cds_features.append(((landmark,cds_id,attributes[6]), (start, end, strand))) # define header for each CDS here
    return cds_features

def load_landmark(fasta_path, landmark):
    genome = Fasta(fasta_path)
    lmk = genome[landmark]
    return lmk


def cds_from_landmark_range(range_start, range_end, cds_list):
    cds_to_return = []
    for tpl in cds_list:
        cds_start, cds_end, _ = tpl[1]
        if (range_start <= cds_start <= range_end and range_start <= cds_end <= range_end ):
            cds_to_return.append(tpl)

    return cds_to_return


def extract_cds_sequences(pyfaidx_landmark, cds_list):
    cds_sequences = []
    for description, location in cds_list:
        start, end, strand = location
        # if strand == '+':
        #     cds_seq = pyfaidx_landmark[start-1:end].seq
        # else:
        #     # TODO Thick about this
        #     cds_seq = pyfaidx_landmark[start-1:end].reverse.seq
        cds_seq = pyfaidx_landmark[start - 1:end].seq
        cds_record = SeqRecord(Seq(cds_seq), id=f'landmark:{description[0]}_ID:{description[1]}_start:{start}_end:{end}_strand:({strand})', description=f'{description[2]};from Amel_HAv3.1')
        cds_sequences.append(cds_record)
    return cds_sequences

def write_to_fasta(seqio_record_list,file_name):
    SeqIO.write(seqio_record_list, file_name, 'fasta')

if __name__ == '__main__':
    genome_path ='/home/daniel/blast/test-data/apis-mellifera/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
    gff_path = '/home/daniel/blast/test-data/apis-mellifera/genomic.gff'
    # specify chromosome / LG 11
    landmark_id = 'NC_037648.1'

    landmark_seq = load_landmark(genome_path, landmark_id)
    cds_all = parse_gff(gff_path,landmark_id)

    # set range of interest. a range since we look at the sequence of yellow-e3 to MRJP9
    yellowe3_start = 2259179
    mrjp9_end = 2325880
    flank_size = 50000

    cds_from_range = cds_from_landmark_range(yellowe3_start - flank_size, mrjp9_end + flank_size, cds_all)
    cds_seqs = extract_cds_sequences(landmark_seq, cds_from_range)

    write_to_fasta(cds_seqs,f'AM-full-gff-attributes-header-CDS-flank-{flank_size}.fasta')