#!/bin/bash

input_filtered_fasta_file="$1"
input_metadata_file="$2"
headers_file="$3"
tmp_file="$4"
output_metadata_file="$5"

grep "^>" "$input_filtered_fasta_file" | cut -d ' ' -f1 | cut -d '>' -f2 > "$tmp_file"

head -n 1 "$headers_file" >> "$output_metadata_file"

zcat "$input_metadata_file" | grep -Fwf "$tmp_file" >> "$output_metadata_file"

echo "WA metadata has been filtered and written to $output_metadata_file."
