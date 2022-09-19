version 1.0

import "tasks/ncov_ingest.wdl" as ncov_ingest

workflow GISAID_INGEST {
  input {
    # ncov ingest
    String GISAID_API_ENDPOINT
    String GISAID_USERNAME_AND_PASSWORD

    # Optionals
    File? cache_nextclade_old
    String? filter  # e.g. "region:Africa" passed to tsv-filters

    Int? cpu
    Int? memory       # in GiB
    Int? disk_size
  }

  call ncov_ingest.gisaid_ingest as ingest {
    input:
      GISAID_API_ENDPOINT = GISAID_API_ENDPOINT,
      GISAID_USERNAME_AND_PASSWORD = GISAID_USERNAME_AND_PASSWORD,

      # optionals
      cache_nextclade_old = cache_nextclade_old,      
      filter = filter,

      cpu = cpu,
      memory = memory,
      disk_size = disk_size
  }

  output {
    # ncov-ingest output either gisaid or genbank
    File sequences_fasta = ingest.sequences_fasta
    File metadata_tsv = ingest.metadata_tsv

    File nextclade_tsv = ingest.nextclade_tsv
    String last_run = ingest.last_run
  }
}
