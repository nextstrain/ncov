## Aim

To build a representative 100k dataset which is available for testing / developing builds locally.
This is intended to run weekly via a GitHub action (which triggers jobs to be run on AWS).
It will upload these files:

* `s3://nextstrain-data/files/ncov/open/100k/metadata.tsv.xz`
* `s3://nextstrain-data/files/ncov/open/100k/sequences.fasta.xz`
* `s3://nextstrain-ncov-private/100k/metadata.tsv.xz`
* `s3://nextstrain-ncov-private/100k/sequences.fasta.xz`

While this profile is not recommended to be run locally, you can see what rules would be run via:

```
snakemake --cores 1 --configfile nextstrain_profiles/100k/config-gisaid.yaml -npf upload --dag | dot -Tpdf > dag-100k-gisaid.pdf
snakemake --cores 1 --configfile nextstrain_profiles/100k/config-open.yaml -npf upload --dag | dot -Tpdf > dag-100k-open.pdf
```

To run manually you can trigger the GitHub action (recommended) or run the jobs locally via:
```
nextstrain build --aws-batch --cpus 16 --memory 31GiB --detach . \
  --configfile nextstrain_profiles/100k/config-gisaid.yaml \
  -f upload
nextstrain build --aws-batch --cpus 16 --memory 31GiB --detach . \
  --configfile nextstrain_profiles/100k/config-open.yaml \
  -f upload
```
