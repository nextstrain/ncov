Test script to sanitize metadata.

  $ pushd "$TESTDIR" > /dev/null

Deduplicate metadata by strain name.
Out of 6 records, 2 are duplicates, so only 4 records plus the header should be retained.
Use a small chunk size to ensure sanitizing works over multiple loops through the metadata.

  $ python3 ../scripts/sanitize_metadata.py \
  >  --metadata unsanitized_metadata.tsv \
  >  --metadata-id-columns 'Virus name' strain \
  >  --database-id-columns gisaid_epi_isl genbank_accession \
  >  --metadata-chunk-size 2 \
  >  --output "$TMP/metadata.tsv"
  $ wc -l "$TMP/metadata.tsv"
  \s*5 .* (re)

Confirm that the duplicate record that we retained has the latest GISAID accession.

  $ grep EPI_ISL_2 "$TMP/metadata.tsv" | wc -l
  \s*1 (re)
  $ rm -f "$TMP/metadata.tsv"

Repeat the process with a tarball, to confirm that we can process these files in the same way.

  $ python3 ../scripts/sanitize_metadata.py \
  >  --metadata unsanitized_metadata.tar.gz \
  >  --metadata-id-columns 'Virus name' strain \
  >  --database-id-columns gisaid_epi_isl genbank_accession \
  >  --output "$TMP/metadata.tsv"
  Extracted metadata file from unsanitized_metadata.tar.gz to .*/unsanitized_metadata.tsv (re)
  Cleaning up temporary files in .* (re)
  $ wc -l "$TMP/metadata.tsv"
  \s*3 .* (re)
  $ rm -f "$TMP/metadata.tsv"

Throw an error when metadata doesn't have any of the requested strain id columns.

  $ python3 ../scripts/sanitize_metadata.py \
  >  --metadata unsanitized_metadata.tsv \
  >  --metadata-id-columns strain \
  >  --database-id-columns gisaid_epi_isl genbank_accession \
  >  --output "$TMP/metadata.tsv"
  ERROR: None of the possible id columns (['strain']) were found in the metadata's columns ('Virus name', 'gender', 'date', 'gisaid_epi_isl')
  [1]

Print a warning when metadata doesn't have any of the requested database id columns.
This should produce a warning about skipping duplicate resolution, an error about the duplicates found, and a file of duplicate strain names.

  $ python3 ../scripts/sanitize_metadata.py \
  >  --metadata unsanitized_metadata.tsv \
  >  --metadata-id-columns "Virus name" \
  >  --database-id-columns genbank_accession \
  >  --output "$TMP/metadata.tsv"
  WARNING: Skipping deduplication of metadata records. None of the possible database id columns (['genbank_accession']) were found in the metadata's columns ('Virus name', 'gender', 'date', 'gisaid_epi_isl')
  ERROR: 2 strains have duplicate records. See 'unsanitized_metadata.tsv.duplicates.txt' for more details.
  [1]
  $ sort unsanitized_metadata.tsv.duplicates.txt
  hCoV-19/LocalVirus/2/2021
  hCoV-19/OneVirus/1/2020
  $ rm -f unsanitized_metadata.tsv.duplicates.txt

Throw an error when metadata contain duplicates.

  $ python3 ../scripts/sanitize_metadata.py \
  >  --metadata unsanitized_metadata.tsv \
  >  --metadata-id-columns 'Virus name' strain \
  >  --database-id-columns gisaid_epi_isl genbank_accession \
  >  --error-on-duplicate-strains \
  >  --output "$TMP/metadata.tsv"
  ERROR: 2 strains have duplicate records. See 'unsanitized_metadata.tsv.duplicates.txt' for more details.
  [1]
  $ sort unsanitized_metadata.tsv.duplicates.txt
  hCoV-19/LocalVirus/2/2021
  hCoV-19/OneVirus/1/2020
  $ rm -f unsanitized_metadata.tsv.duplicates.txt

Rename fields and strip prefixes.

  $ head -n 1 "unsanitized_metadata.tsv"
  Virus name\tgender\tdate\tgisaid_epi_isl (esc)
  $ grep "hCoV-19" "unsanitized_metadata.tsv" | wc -l
  \s*5 (re)
  $ grep "SARS-CoV-2" "unsanitized_metadata.tsv" | wc -l
  \s*1 (re)
  $ python3 ../scripts/sanitize_metadata.py \
  >  --metadata unsanitized_metadata.tsv \
  >  --rename-fields "Virus name=strain" "gender=sex" \
  >  --strip-prefixes "hCoV-19/" "SARS-CoV-2/" \
  >  --output "$TMP/metadata.tsv"
  $ head -n 1 "$TMP/metadata.tsv"
  strain\tsex\tdate\tgisaid_epi_isl (esc)
  $ grep "hCoV-19" "$TMP/metadata.tsv"
  [1]
  $ grep "SARS-CoV-2" "$TMP/metadata.tsv"
  [1]
  $ rm -f "$TMP/metadata.tsv"
