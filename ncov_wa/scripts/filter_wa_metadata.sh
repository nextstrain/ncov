#!/bin/bash
input_metadata_file="$1"
output_metadata_file="$2"

# Extract header line and write to output
zcat "$input_metadata_file" | head -n 1 > "$output_metadata_file"

# Filter rows and append to output
zcat "$input_metadata_file" | \
awk -F'\t' -v strain_idx="$strain_idx" -v division_idx="$division_idx" 'NR > 1 && ($strain_idx ~ /USA\/WA/ || $division_idx == "Washington" || $division_idx == "WA")' \
>> "$output_metadata_file"
echo "Filtered metadata written to $output_metadata_file"
