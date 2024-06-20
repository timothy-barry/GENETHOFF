from Bio import SeqIO
import gzip
import sys

def merge_fastq_files(input_file1, input_file2, output_file):
    with gzip.open(input_file1, "rt") as f1, gzip.open(input_file2, "rt") as f2, gzip.open(output_file, "wt") as out_f:
        for record1, record2 in zip(SeqIO.parse(f1, "fastq"), SeqIO.parse(f2, "fastq")):
            merged_sequence = record1.seq + record2.seq
            merged_quality = record1.letter_annotations["phred_quality"] + record2.letter_annotations["phred_quality"]

            merged_record = record1
            del merged_record.letter_annotations["phred_quality"]
            merged_record.seq = merged_sequence
            merged_record.letter_annotations["phred_quality"] = merged_quality

            SeqIO.write(merged_record, out_f, "fastq")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_fastq.py input_file1 input_file2 output_file")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]
    output_file = sys.argv[3]

    merge_fastq_files(input_file1, input_file2, output_file)
    print("FASTQ files merged successfully.")

