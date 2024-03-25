from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# same for all blast types e.g. blastn, tblastx etc; note the -outfmt 6 must be used
def read_blast_results_to_list(file_path):
    content = []
    with open(file_path, 'r') as f:
        for line in f:
           content.append(tuple(line.split()))

    return content

def get_unique_landmarks(content_touple_list):
    landmarks = [t[0] for t in content_touple_list]
    return list(set(landmarks))

def get_landmark_SeqIO_records(pyfaidx_genome, list_of_landmarks):
    record_list = []
    for lmk in list_of_landmarks:
        lmk_seq = pyfaidx_genome[lmk][:].seq
        lmk_header = pyfaidx_genome[lmk].long_name
        record = SeqRecord(Seq(lmk_seq),id=lmk_header)
        record_list.append(record)

    return record_list

# TODO genome id needed to identify origin
def get_tblastx_match_SeqIO_records(pyfaidx_genome,tblastx_content):
    records = []
    current_seq = ''
    file_name = pyfaidx_genome.filename.split('/')[-1]
    genome_version = '_'.join(file_name.split('_')[2:4])
    for match in tblastx_content:
        landmark = match[0]
        start = int(match[6])
        end = int(match[7])
        # pyfaidx slicing: start is not included but end is -> need to subract 1 from start to include
        if end < start:
            current_seq = pyfaidx_genome[landmark][end-1:start].reverse.seq
            header = f'>{landmark}_start:{end}_end:{start}_(-)'
            record = SeqRecord(Seq(current_seq),id=header,description=f'Genome: {genome_version}')
            records.append(record)
        else:
            current_seq = pyfaidx_genome[landmark][start-1:end].seq
            header = f'>{landmark}_start:{start}_end:{end}_(+)'
            record = SeqRecord(Seq(current_seq),id=header,description=f'Genome: {genome_version}')
            records.append(record)

    return records



def write_to_fasta(seqio_record_list, file_name):
    SeqIO.write(seqio_record_list, file_name, 'fasta')

if __name__ == '__main__':
    file_path = '/home/daniel/blast/res-af'
    blastn_content = read_blast_results_to_list(file_path)
    unique_landmarks = get_unique_landmarks(blastn_content)
    genome_path = '/home/daniel/blast/test-data/apis-florea/GCF_000184785.3_Aflo_1.1_genomic.fna'
    genome = Fasta(genome_path)
    records = get_landmark_SeqIO_records(genome, unique_landmarks)
    print(len(records[0].seq))
    write_to_fasta(records, 'apis-florea_ROI.fasta')

    tblastx_results = '/home/daniel/blast/resx-af'
    tblastx_content = read_blast_results_to_list(tblastx_results)
    tblastx_records = get_tblastx_match_SeqIO_records(genome,tblastx_content)
    write_to_fasta(tblastx_records,'apis-florea_tblastx-selection.fasta')

