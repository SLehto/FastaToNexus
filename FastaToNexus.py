import os
import sys
from Bio import SeqIO

def read_fasta_files(fasta_files):
    taxa = set()
    alignments = {}
    
    for fasta_file in fasta_files:
        alignment = {}
        seq_lengths = set()
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            taxa.add(record.id)
            seq_length = len(record.seq)
            seq_lengths.add(seq_length)
            alignment[record.id] = str(record.seq)
        
        if len(seq_lengths) > 1:
            print(f"Error: Sequences in {fasta_file} have different lengths!")
            sys.exit(1)
        
        alignments[fasta_file] = alignment
    
    return taxa, alignments

def get_sequence_lengths(alignments):
    return {file: len(next(iter(align.values()), '')) for file, align in alignments.items()}

def pad_missing_sequences(taxa, alignments, seq_lengths):
    concatenated_sequences = {taxon: "" for taxon in taxa}
    partitions = {}
    start_pos = 1
    
    for file, alignment in alignments.items():
        length = seq_lengths[file]
        end_pos = start_pos + length - 1
        partitions[file] = (start_pos, end_pos)
        
        for taxon in taxa:
            concatenated_sequences[taxon] += alignment.get(taxon, "N" * length)
        
        start_pos = end_pos + 1
    
    return concatenated_sequences, partitions

def write_nexus(output_file, taxa, concatenated_sequences, partitions):
    with open(output_file, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={len(taxa)} NCHAR={len(next(iter(concatenated_sequences.values())))};\n")
        f.write("  FORMAT DATATYPE=DNA MISSING=N GAP=-;\n")
        f.write("  MATRIX\n")
        
        for taxon in sorted(taxa):
            f.write(f"  {taxon}  {concatenated_sequences[taxon]}\n")
        
        f.write("  ;\n")
        f.write("END;\n\n")
        
        f.write("BEGIN SETS;\n")
        for file, (start, end) in partitions.items():
            partition_name = os.path.basename(file).replace(".", "_")
            f.write(f"  CHARSET {partition_name} = {start}-{end};\n")
        f.write("END;\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py output.nex input1.fasta input2.fasta ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    fasta_files = sys.argv[2:]
    
    taxa, alignments = read_fasta_files(fasta_files)
    seq_lengths = get_sequence_lengths(alignments)
    concatenated_sequences, partitions = pad_missing_sequences(taxa, alignments, seq_lengths)
    write_nexus(output_file, taxa, concatenated_sequences, partitions)
    
    print(f"NEXUS file '{output_file}' successfully created.")

if __name__ == "__main__":
    main()