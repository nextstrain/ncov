#!/bin/bash
input_fasta_file="$1" #Nextstrain remote dataset file
filtered_metadata_file="$2" #WA metadata file
tmp_file="$3"
output_fasta_file="$4" #Filtered WA fasta file

# Decompress fasta files into temp
xzcat "$input_fasta_file" > temp_decompressed.fasta

# 1. Extract strain names (assumes 'strain' is the first column)
cut -f1 "$filtered_metadata_file" | tail -n +2 | tr -d '\r' > "$tmp_file" #ensures temp files does not contain and \r characters
tmp_file=$(echo "$tmp_file" | tr -d '\r')

#matches strain name to fasta header
seqtk subseq temp_decompressed.fasta "$tmp_file" > "$output_fasta_file"
echo "WA seqs have been filtered and written to $output_fasta_file."
