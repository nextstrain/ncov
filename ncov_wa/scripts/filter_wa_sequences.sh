#!/bin/bash

input_file="$1"
output_file="$2"
pattern="USA/WA"

xzcat "$input_file" | awk -v pattern="$pattern" '
BEGIN { RS=">"; ORS="" }
NR > 1 && $0 ~ pattern {
  split($0, header_seq, "\n");
  header = header_seq[1];
  seq = "";
  for (i = 2; i <= length(header_seq); i++) {
    seq = seq header_seq[i];
  }
  if (length(seq) >= 27000) {
    print ">" header "\n" seq "\n";
  }
}' > "$output_file"


echo "WA seqs have been filtered and written to $output_file."
