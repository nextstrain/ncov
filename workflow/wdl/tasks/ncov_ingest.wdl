version 1.0

task gisaid_ingest {
  input {
    String GISAID_API_ENDPOINT
    String GISAID_USERNAME_AND_PASSWORD

    # String? AWS_DEFAULT_REGION
    # String? AWS_ACCESS_KEY_ID
    # String? AWS_SECRET_ACCESS_KEY
    # String? SLACK_TOKEN
    # String? SLACK_CHANNEL

    # Optionals
    File? cache_nextclade_old
    String? filter  # e.g. "region:Africa" passed to tsv-filters

    String giturl = "https://github.com/nextstrain/ncov-ingest/archive/refs/heads/master.zip"

    Int cpu = if ( defined(cache_nextclade_old) ) then 16 else 96
    Int disk_size = 1500  # In GiB
    Float memory = if ( defined(cache_nextclade_old) ) then 64 else 180
  }

  command <<<
    echo "cpu: ~{cpu}; disk_size: ~{disk_size}; memory: ~{memory}"
    nextclade --version
    zstd --version

    # Set up env variables
    export GISAID_API_ENDPOINT="~{GISAID_API_ENDPOINT}"
    export GISAID_USERNAME_AND_PASSWORD="~{GISAID_USERNAME_AND_PASSWORD}"

    export PROC=`nproc`
    export temp_mem="~{memory}"
    export MEM=${temp_mem%.*}000

    # Pull ncov-ingest repo
    wget -O master.zip ~{giturl}
    NCOV_INGEST_DIR=`unzip -Z1 master.zip | head -n1 | sed 's:/::g'`
    unzip master.zip

    # Try to use cache to shorten runtime
    if [[ -n "~{cache_nextclade_old}" ]]; then
      export NEXTCLADE_CACHE="~{cache_nextclade_old}"

      # Detect and decompress xz or zst files
      if [[ $NEXTCLADE_CACHE == *.xz ]]; then
        cp ~{cache_nextclade_old} .
        xz -T0 --decompress ~{basename(select_first([cache_nextclade_old,'']))}
        export NEXTCLADE_CACHE="~{basename(select_first([cache_nextclade_old,'']),'.xz')}"
      elif [[ $NEXTCLADE_CACHE == *.zst ]]; then
        cp ~{cache_nextclade_old} .
        zstd -T0 -d ~{basename(select_first([cache_nextclade_old,'']))}
        export NEXTCLADE_CACHE="~{basename(select_first([cache_nextclade_old,'']),'.zst')}"
      fi

      cp $NEXTCLADE_CACHE ${NCOV_INGEST_DIR}/data/gisaid/nextclade_old.tsv
      echo "nextclade_old.tsv has " `wc -l ${NCOV_INGEST_DIR}/data/gisaid/nextclade_old.tsv` " records."
    fi

    # Navigate to ncov-ingest directory, and call snakemake
    cd ${NCOV_INGEST_DIR}

    nextstrain build \
      --native \
      --cpus $PROC \
      --memory ~{memory}GiB \
      --exec env \
      . \
        snakemake \
          --configfile config/local_gisaid.yaml \
          --cores $PROC \
          --resources mem_mb=$MEM \
          --printshellcmds

    # === prepare output
    # Date stamp last run (YYYY-MM-DD)
    LAST_RUN=`date +%F`

    cd ..
    echo ${LAST_RUN} > LAST_RUN
    ls -l ${NCOV_INGEST_DIR}/data/*
    echo "Number of sequences: " `wc -l ${NCOV_INGEST_DIR}/data/gisaid/metadata_transformed.tsv`

    # Optional filter
    if [[ -n "~{filter}" ]]; then
      cat ${NCOV_INGEST_DIR}/data/gisaid/metadata_transformed.tsv \
      | tsv-filter -H --str-in-fld ~{filter} \
      | zstd -T0 -o gisaid_metadata.tsv.zst

      echo "With ~{filter}, filtered down to: " `zstd -d -c gisaid_metadata.tsv.zst | wc -l`

      zstd -T0 -d -c gisaid_metadata.tsv.zst \
      | tsv-select -H -f 'virus','date','date_submitted' \
      | sed 1d \
      | awk -F "\t" '{ print $1"|"$2"|"$3 }' > strain_list.txt

      wget https://raw.githubusercontent.com/santiagosnchez/faSomeRecords/master/faSomeRecords.py

      cat ${NCOV_INGEST_DIR}/data/gisaid/sequences.fasta \
      | python faSomeRecords.py --fasta /dev/stdin --list strain_list.txt --stdout \
      | zstd -T0 -o gisaid_sequences.fasta.zst

    else
      mv ${NCOV_INGEST_DIR}/data/gisaid/sequences.fasta gisaid_sequences.fasta
      mv ${NCOV_INGEST_DIR}/data/gisaid/metadata_transformed.tsv gisaid_metadata.tsv
      zstd -T0 gisaid_sequences.fasta
      zstd -T0 gisaid_metadata.tsv
    fi

    # prepare output caches
    if [[ -f "${NCOV_INGEST_DIR}/data/gisaid/nextclade.tsv" ]]; then
      mv ${NCOV_INGEST_DIR}/data/gisaid/nextclade.tsv gisaid_nextclade.tsv
    elif [[ -f ${NCOV_INGEST_DIR}/data/gisaid/nextclade_old.tsv ]]; then
      mv ${NCOV_INGEST_DIR}/data/gisaid/nextclade_old.tsv gisaid_nextclade.tsv
    else
      touch gisaid_nextclade.tsv
    fi
    zstd -T0 gisaid_nextclade.tsv

    ls -lh gisaid_nextclade.tsv.zst gisaid_sequences.fasta.zst gisaid_metadata.tsv.zst
  >>>

  output {
    # Ingested gisaid sequence and metadata files
    File sequences_fasta = "gisaid_sequences.fasta.zst"
    File metadata_tsv = "gisaid_metadata.tsv.zst"

    # cache for next run
    File nextclade_tsv = "gisaid_nextclade.tsv.zst"
    String last_run = read_string("LAST_RUN")
  }

  runtime {
    docker: "nextstrain/ncov-ingest:latest"
    cpu : cpu
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

}

task genbank_ingest {
  input {
    # Optionals
    File? cache_nextclade_old
    String? filter  # e.g. "region:Africa" passed to tsv-filters

    String giturl = "https://github.com/nextstrain/ncov-ingest/archive/refs/heads/master.zip"

    Int cpu = if ( defined(cache_nextclade_old) ) then 16 else 96
    Int disk_size = 1500  # In GiB
    Float memory = if ( defined(cache_nextclade_old) ) then 64 else 180
  }

  command <<<
    echo "cpu: ~{cpu}; disk_size: ~{disk_size}; memory: ~{memory}"
    nextclade --version
    zstd --version

    # Set up env variables
    export PROC=`nproc`
    export temp_mem="~{memory}"
    export MEM=${temp_mem%.*}000

    # Pull ncov-ingest repo
    wget -O master.zip ~{giturl}
    NCOV_INGEST_DIR=`unzip -Z1 master.zip | head -n1 | sed 's:/::g'`
    unzip master.zip

    # Try to use cache to shorten runtime
    if [[ -n "~{cache_nextclade_old}" ]]; then
      export NEXTCLADE_CACHE="~{cache_nextclade_old}"

      # Detect and decompress xz or zst files
      if [[ $NEXTCLADE_CACHE == *.xz ]]; then
        cp ~{cache_nextclade_old} .
        xz -T0 --decompress ~{basename(select_first([cache_nextclade_old,'']))}
        export NEXTCLADE_CACHE="~{basename(select_first([cache_nextclade_old,'']),'.xz')}"
      elif [[ $NEXTCLADE_CACHE == *.zst ]]; then
        cp ~{cache_nextclade_old} .
        zstd -T0 -d ~{basename(select_first([cache_nextclade_old,'']))}
        export NEXTCLADE_CACHE="~{basename(select_first([cache_nextclade_old,'']),'.zst')}"
      fi

      cp $NEXTCLADE_CACHE ${NCOV_INGEST_DIR}/data/genbank/nextclade_old.tsv
      echo "nextclade_old.tsv has " `wc -l ${NCOV_INGEST_DIR}/data/genbank/nextclade_old.tsv` " records."
    fi

    # Navigate to ncov-ingest directory, and call snakemake
    cd ${NCOV_INGEST_DIR}

    nextstrain build \
      --native \
      --cpus $PROC \
      --memory ~{memory}GiB \
      --exec env \
      . \
        snakemake \
          --configfile config/local_genbank.yaml \
          --cores $PROC \
          --resources mem_mb=$MEM \
          --printshellcmds

    # === prepare output
    # Date stamp last run (YYYY-MM-DD)
    LAST_RUN=`date +%F`

    cd ..
    echo ${LAST_RUN} > LAST_RUN
    ls -l ${NCOV_INGEST_DIR}/data/*
    echo "Number of sequences: " `wc -l ${NCOV_INGEST_DIR}/data/genbank/metadata_transformed.tsv`

    # Optional filter
    if [[ -n "~{filter}" ]]; then
      cat ${NCOV_INGEST_DIR}/data/genbank/metadata_transformed.tsv \
      | tsv-filter -H --str-in-fld ~{filter} \
      | zstd -T0 -o genbank_metadata.tsv.zst

      echo "With ~{filter}, filtered down to: " `zstd -d -c genbank_metadata.tsv.zst | wc -l`

      zstd -T0 -d -c genbank_metadata.tsv.zst \
      | tsv-select -H -f 'virus','date','date_submitted' \
      | sed 1d \
      | awk -F "\t" '{ print $1"|"$2"|"$3 }' > strain_list.txt

      wget https://raw.githubusercontent.com/santiagosnchez/faSomeRecords/master/faSomeRecords.py

      cat ${NCOV_INGEST_DIR}/data/genbank/sequences.fasta \
      | python faSomeRecords.py --fasta /dev/stdin --list strain_list.txt --stdout \
      | zstd -T0 -o genbank_sequences.fasta.zst

    else
      mv ${NCOV_INGEST_DIR}/data/genbank/sequences.fasta genbank_sequences.fasta
      mv ${NCOV_INGEST_DIR}/data/genbank/metadata_transformed.tsv genbank_metadata.tsv
      zstd -T0 genbank_sequences.fasta
      zstd -T0 genbank_metadata.tsv
    fi

    # prepare output caches
    if [[ -f "${NCOV_INGEST_DIR}/data/genbank/nextclade.tsv" ]]; then
      mv ${NCOV_INGEST_DIR}/data/genbank/nextclade.tsv genbank_nextclade.tsv
    elif [[ -f ${NCOV_INGEST_DIR}/data/gisaid/nextclade_old.tsv ]]; then
      mv ${NCOV_INGEST_DIR}/data/gisaid/nextclade_old.tsv genbank_nextclade.tsv
    else
      touch genbank_nextclade.tsv
    fi
    zstd -T0 genbank_nextclade.tsv

    ls -lh genbank_nextclade.tsv.zst genbank_sequences.fasta.zst genbank_metadata.tsv.zst
  >>>

  output {
    # Ingested genbank sequence and metadata files
    File sequences_fasta = "genbank_sequences.fasta.zst"
    File metadata_tsv = "genbank_metadata.tsv.zst"

    # cache for next run
    File nextclade_tsv = "genbank_nextclade.tsv.zst"
    String last_run = read_string("LAST_RUN")
  }

  runtime {
    docker: "nextstrain/ncov-ingest:latest"
    cpu : cpu
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

}