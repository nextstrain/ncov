## Omicron only GISAID build

### Prerequisites:

Appropriate AWS credentials

### How to run

**locally**
```
snakemake --cores 4 --configfile my_profiles/omicron_gisaid/config.yaml -p -f auspice/ncov_omicron.json
```

**AWS**

```
nextstrain build --aws-batch --cpus 16 --memory 14GiB --detach . --configfile my_profiles/omicron_gisaid/config.yaml -p -f auspice/ncov_omicron.json
```