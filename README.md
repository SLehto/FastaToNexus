This is a python script to concatenate multiple (pre-aligned) FASTA files into a NEXUS-formatted file. It reads multiple FASTA files, identifies unique taxa, verifies that each sequence in FASTA are equally long, concatenates the sequences, and outputs a NEXUS-formatted file with a sets block indicating the partitions corresponding to the original FASTA files.

1. Reads multiple FASTA files and extracts unique taxa and sequences.
2. Checks that all sequences within each FASTA file have the same length.
3. Identifies the length of each alignment and ensures missing sequences are filled with gaps (N).
4. Concatenates the sequences and tracks partitions for the sets block.
5. Outputs a NEXUS file.

Usage:
python FastaToNexus.py output.nex file1.fasta file2.fasta file3.fasta

This script ensures that taxa missing from some alignments are properly padded while preserving partitions for downstream phylogenetic analysis.
 
