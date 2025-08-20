#!/bin/bash

datasets download virus genome taxon sars-cov-2 --usa-state WA --filename wa_sars2.zip

unzip wa_sars2.zip -d wa_sars2

cp wa_sars2/ncbi_dataset/data/genomic.fna wa_sars2/ncbi_dataset/data/wa_sequences.fasta

dataformat tsv virus-genome --inputfile wa_sars2/ncbi_dataset/data/data_report.jsonl --fields accession,sourcedb,isolate-lineage,geo-region,geo-location,isolate-collection-date,release-date,update-date,length,host-name,is-lab-host,isolate-lineage-source,bioprojects,biosample-acc,sra-accs,submitter-names,submitter-affiliation > wa_sars2/ncbi_dataset/data/wa_metadata.tsv

rscript NCBI_to_GISAID_format.R

echo "NCBI WA metadata and sequences downloaded and cleaned"
