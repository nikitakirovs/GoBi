import os
import argparse
import fastaparser

# Input FASTA files directory
input_dir = "input_fasta_files"
# Output directory for the results
output_dir = "blastn_results"
# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Specify the path to your FASTA file
fasta_file_path = "data/sequences.fna"

#creat new sequence files for blast
file1 = open("input_fasta_files/seq1.txt", "w")

# Read the FASTA file
with open(fasta_file_path) as fasta_file:
    parser = fastaparser.Reader(fasta_file)

    # create the db for blastn
    db = "data/Bombus_impatiens.fa"
    cmd1 = "makeblastdb -in " + db + " -dbtype nucl"
    os.system(cmd1)


    # Iterate through the sequences
    for seq in parser:
        print(f"ID: {seq.id}")
        print(f"Description: {seq.description}")
        print(seq.sequence_as_string())
        file1.write(seq.sequence_as_string())
        filepath = fasta_file_path
        output = output_dir + "/" + "output.txt"
        cmd2 = "blastn -query " + filepath + " -db " + db + " -out " + output + " -outfmt 7"
        os.system(cmd2)


file1.close()

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("--db")
parser.add_argument("--file")
args = parser.parse_args()


# or simply blast the whole fasta file:
# blastn -query data/seq1.fna -db data/Apis_florea.fa -out output.txt
