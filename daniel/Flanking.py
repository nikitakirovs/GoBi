from Bio import SeqIO
import re



def extract_id_location_raw_to_dict(seqrecord_obj):
    all_descriptions_str_list = list()
    for record in seqrecord_obj:
        all_descriptions_str_list.append(record.description)

    id_location_dict = {}

    for description in all_descriptions_str_list:
        description_split = description.split()
        id_location_dict[description_split[0]] = description_split[-2]

    return id_location_dict

# takes id_location_dict returned from extract_id_location_to_dict
def location_contains_order_operator(id_location_raw_dict):
    contains_order = False
    for location in id_location_raw_dict.values():
        if 'order' in location:
            contains_order = True
        if contains_order:
            break
    return contains_order

"""
filter out all CDS where location contains 'order' operator, aka iffy ones:

order(location,location, ... location) 
The elements can be found in the 
specified order (5' to 3' direction), but nothing is implied about the 
reasonableness about joining them

or CDS which are on the complementary DNA strand
"""

def extract_uninterrupted_and_join_cds_locations_to_dict(id_location_dict_raw):
    pattern = re.compile('\[location=[\.0-9]+\]|\[location=join\([\.,0-9]+\)\]')
    location_dict_processed = {}
    for i,l in id_location_dict_raw.items():
        if pattern.match(l):
            location_dict_processed[i] = extract_start_and_end_index(l)
    return location_dict_processed


# extract start and end index. can be simple range e.g 134..235 or join(range1,range2,...,range n)
def extract_start_and_end_index(location_string):
    pattern = re.compile('([0-9]+)\.\.([0-9]+|\,|\.\.)+')
    matches = pattern.findall(location_string) # list of matches, tested that it always contains only one match because of regex structure
    strx, stry = matches[0]
    return (int(strx), int(stry))


def get_location_of_gene(some_gene_id_string,location_dict_processed):
    return location_dict_processed[some_gene_id_string]

# takes location of one CDS or some other index range -> range_tuple
def in_range(cds_location_to_check, search_range):
    x,y = cds_location_to_check
    # we can assume x < y
    if search_range[0] <= x <=search_range[-1] and search_range[0] <= y <= search_range[-1]:
        return True
    else:
        return False

def get_flanking_genes_of_one_cds_symmetric_offset(some_gene_id_string, offset_int, location_dict_processed):
    gene_start, gene_end = get_location_of_gene(some_gene_id_string,location_dict_processed)
    flank_window_left = (gene_start - offset_int, gene_start)
    flank_window_right = (gene_end, gene_end + offset_int)

    flanking_left = {}
    flanking_right = {}
    
    for i,l in location_dict_processed.items():
        if in_range(l,flank_window_left):
            flanking_left[i] = l
        elif in_range(l,flank_window_right):
            flanking_right[i] = l
        else:
            continue

    return [flanking_left,flanking_right]


def get_cds_in_range(range_tuple, location_dict_processed):
    left, right = range_tuple
    in_range_dict = {}
    for i, l in location_dict_processed.items():
        if in_range(l, range_tuple):
            in_range_dict[i] = l
        else:
            continue

    return in_range_dict

def get_flanking_genes_of_specific_range_symmetric_offset(range_tuple_int, offset_int, location_dict_processed):
    left, right = range_tuple_int
    flank_window_left = (left - offset_int, left)
    flank_window_right = (right, right + offset_int)

    flanking_left = {}
    flanking_right = {}

    for i, l in location_dict_processed.items():
        if in_range(l, flank_window_left):
            flanking_left[i] = l
        elif in_range(l, flank_window_right):
            flanking_right[i] = l
        else:
            continue

    return [flanking_left, flanking_right]


def extract_flanking_records_to_list(seqrecord_obj, list_flanking_left_flanking_right):
    flanking_records = []
    flanking_left, flanking_right = list_flanking_left_flanking_right
    id_list = list(flanking_left.keys()) + list(flanking_right.keys())
    for record in seqrecord_obj:
        if record.id in id_list:
            flanking_records.append(record)

    return flanking_records

def extract_cds_in_range_records(seqrecord_obj, in_range_dict):
    ids_only = list(in_range_dict.keys())
    cds_in_range_records = []
    for record in seqrecord_obj:
        if record.id in ids_only:
            cds_in_range_records.append(record)

    return cds_in_range_records

def write_records_to_fasta(flanking_records, file_name):
    SeqIO.write(flanking_records, file_name, 'fasta')

def remove_records_with_seqs_shorter_than(length_int, record_list):
    result_list = []
    for record in record_list:
        if len(record.seq) < length_int:
            continue
        else:
            result_list.append(record)

    return result_list

# if this file is not imported as a module run the following code
if __name__ == '__main__':
    file_path = "/home/me/gobi/apis-mellifera-cds/data/GCF_003254395.2/cds_from_genomic.fna"
    seq_record = SeqIO.parse(file_path, "fasta")
    id_loc_dict_raw = extract_id_location_raw_to_dict(seq_record)
    id_loc_dict_processed = extract_uninterrupted_and_join_cds_locations_to_dict(id_loc_dict_raw)
    some_id = 'lcl|NC_037648.1_cds_NP_001091698.1_16484'
    flanking_left, flanking_right = get_flanking_genes_of_one_cds_symmetric_offset(some_id,50000,id_loc_dict_processed)
    # print(f'Gene: {some_id}')
    # print(f'Flanking CDS TOTAL Count on (+) strand: {len(flanking_left)+len(flanking_right)}')
    # print(f'Flanking LEFT Count: {len(flanking_left)}')
    # print(f'Flanking RIGHT Count: {len(flanking_right)}')
    # print(f'Flanking CDS on (+) strand LEFT:')
    # for i in flanking_left:
    #     print(f'\t{i}')
    #
    # print(f'\nFlanking CDS on (+) strand RIGHT:\n')
    # for i in flanking_right:
    #     print(f'\t{i}')
    filter_length = 500
    print(f'filter length: {filter_length}')
    flanking_records = extract_flanking_records_to_list(SeqIO.parse(file_path,'fasta'),[flanking_left,flanking_right])
    flanking_records = remove_records_with_seqs_shorter_than(filter_length, flanking_records)
    print(f'flanking remaining: {len(flanking_records)}')
    write_records_to_fasta(flanking_records, 'flanking.fasta')

    cds_in_range_dict = get_cds_in_range((2259179,2325880),id_loc_dict_processed)
    cds_in_range_records = extract_cds_in_range_records(SeqIO.parse(file_path,'fasta'), cds_in_range_dict)
    cds_in_range_records = remove_records_with_seqs_shorter_than(filter_length,cds_in_range_records)
    print(f'main genes remaining: {len(cds_in_range_records)}')
    write_records_to_fasta(cds_in_range_records,'genes.fasta')