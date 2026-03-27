Test script to annotate metadata with sequence index.

  $ pushd "$TESTDIR" > /dev/null

Annotate metadata with an index where strains have numeric names.

  $ python3 ../scripts/annotate_metadata_with_index.py \
  >  --metadata numeric_strains/metadata.tsv \
  >  --sequence-index numeric_strains/sequence_index.tsv \
  >  --output "$TMP/metadata_with_index.tsv"
  $ wc -l "$TMP/metadata_with_index.tsv"
  \s*3 .* (re)

  $ rm -f "$TMP/metadata_with_index.tsv"

Annotate metadata with an index where metadata uses the "name" column for strain names.
This should produce the same result as above, except with a "name" column in the output.

  $ python3 ../scripts/annotate_metadata_with_index.py \
  >  --metadata numeric_strains/metadata_by_name.tsv \
  >  --sequence-index numeric_strains/sequence_index.tsv \
  >  --output "$TMP/metadata_with_index.tsv"
  $ wc -l "$TMP/metadata_with_index.tsv"
  \s*3 .* (re)

  $ rm -f "$TMP/metadata_with_index.tsv"

  $ popd > /dev/null
