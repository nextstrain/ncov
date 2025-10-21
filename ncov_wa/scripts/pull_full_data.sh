#!/bin/bash

full_metadata_url="$1"
full_sequences_url="$2"
full_metadata="$3"
full_sequences="$4"

wget -O "$full_metadata" "$full_metadata_url"
chmod +x "$full_metadata"

wget -O "$full_sequences" "$full_sequences_url"
chmod +x "$full_sequences"

echo "Full SARS-CoV-2 remote dataset from Nextstrain is complete."
