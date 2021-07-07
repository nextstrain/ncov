Integration tests for nCoV pipeline.

Note that running these tests requires setup steps, and that each test can only
run one-at-a-time due to the shared use of the test environment as otherwise
snakemake may use intermediate files from previous runs, thus producing
inconsistent test results.

Cram should be run in an environment which can run the pipeline via
`cram --preserve-env tests/different-inputs.t` or similar.

Set-up test environment. We could set up the correct data inside $TMP for each test
if we prefer. For simplicity, we create a directory "output".

  $ pushd "$TESTDIR" > /dev/null
  $ basename $( pwd )
  tests
  $ rm -rf output && mkdir output && cd output
  $ cp -r ../../defaults . && cp -r ../../scripts . && mkdir data/ && cp ../../data/references* data/
  $ cd ../..
  $ basename $( pwd )
  ncov

Test various input starting points, all from local (.xz) compressed files

  $ snakemake --directory tests/output --profile tests/local-inputs-compressed \
  > auspice/ncov_test-local-compressed.json >tests/output/local-inputs-compressed.cram.log.txt 2>&1

  $ python3 tests/check_auspice_json.py --json tests/output/auspice/ncov_test-local-compressed.json \
  > --attr region --values "North America" "Europe" "Asia" "Oceania"

  $ rm -rf tests/output/results

Test various input starting points, all from remote (.xz) compressed files

  $ snakemake --directory tests/output --profile tests/remote-inputs-compressed \
  > auspice/ncov_test-remote-compressed.json >tests/output/remote-inputs-compressed.cram.log.txt 2>&1

  $ python3 tests/check_auspice_json.py --json tests/output/auspice/ncov_test-remote-compressed.json \
  > --attr region --values "North America" "Europe" "Asia" "Oceania"

  $ rm -rf tests/output/results tests/output/data/downloaded_test*compressed*

Test various input starting points, all from local uncompressed files

  $ cp tests/local-inputs-compressed/data/*xz tests/local-inputs-uncompressed/data/

  $ for i in tests/local-inputs-uncompressed/data/*.xz; do xz -d $i; done

  $ snakemake --directory tests/output --profile tests/local-inputs-uncompressed \
  > auspice/ncov_test-local-uncompressed.json >tests/output/local-inputs-uncompressed.cram.log.txt 2>&1

  $ python3 tests/check_auspice_json.py --json tests/output/auspice/ncov_test-local-uncompressed.json \
  > --attr region --values "North America" "Europe" "Asia" "Oceania"

  $ rm -rf tests/output/results tests/local-inputs-uncompressed/data/*.fasta tests/local-inputs-uncompressed/data/*.tsv

Test various input starting points which support remote uncompressed files (this is a subset of available inputs)

  $ snakemake --directory tests/output  --profile tests/remote-inputs-uncompressed \
  > auspice/ncov_test-remote-uncompressed.json >tests/output/remote-inputs-uncompressed.cram.log.txt 2>&1

  $ python3 tests/check_auspice_json.py --json tests/output/auspice/ncov_test-remote-uncompressed.json \
  > --attr region --values "North America" "Europe" "Asia" "Oceania"

  $ rm -rf tests/output/results data/downloaded_test*uncompressed*