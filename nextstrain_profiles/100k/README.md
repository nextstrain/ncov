## Aim

To build a representative 100k dataset which is available for testing / developing builds locally.
This is intended to run weekly via a GitHub action (which triggers a job to be run on AWS).
It will make two files available:

* `s3://nextstrain-ncov-private/100k/metadata.tsv.xz`
* `s3://nextstrain-ncov-private/100k/sequences.fasta.xz`

While this profile is not recommended to be run locally, you can see what rules would be run via:

```
snakemake --cores 1 --configfile nextstrain_profiles/100k/config.yaml -npf upload --dag | dot -Tpdf > dag.pdf
```

To run manually you can trigger the GitHub action or run the job locally via:
```
nextstrain build --aws-batch --cpus 16 --memory 31GiB --detach . \
  --configfile nextstrain_profiles/100k/config.yaml \
  -f upload
```
