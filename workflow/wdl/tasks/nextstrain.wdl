version 1.0

task nextstrain_build {
  input {
    File? sequence_fasta
    File? metadata_tsv
    File? context_targz  #<= optional contextual sequence
    String build_name = "example"

    File? configfile_yaml # e.g. builds.yaml
    File? custom_zip      # <= since custom is private
    # String? custom_url = "path to public github"  # Our custom config files are private
    String? active_builds # Wisconsin,Minnesota,Washington

    String? s3deploy      # "s3://nextstrain-staging/"
    String? AWS_ACCESS_KEY_ID
    String? AWS_SECRET_ACCESS_KEY

    String pathogen_giturl = "https://github.com/nextstrain/ncov/archive/refs/heads/master.zip"
    Int cpu = 8
    Int disk_size = 30  # In GiB.  Could also check size of sequence or metadata files
    Float memory = 3.5
  }
  command <<<
    # Pull ncov, zika or similar pathogen repo
    wget -O master.zip ~{pathogen_giturl}
    INDIR=`unzip -Z1 master.zip | head -n1 | sed 's:/::g'`
    unzip master.zip

    # If a config file (builds.yaml) file is not provided, generate one
    export CONFIGFILE_FLAG=""
    if [[ -n "~{sequence_fasta}" ]]; then
      if [ -z "~{metadata_tsv}" ]; then
        echo "Error: Provided sequence: ~{sequence_fasta} but missing metadata tsv file."
        exit 1
      fi

      if [[ -z "~{configfile_yaml}" ]]; then
    cat << EOF > builds.yaml
    inputs:
    - name: "~{build_name}"
      metadata: ~{metadata_tsv}
      sequences: ~{sequence_fasta}
    - name: opengenbank
      metadata: https://data.nextstrain.org/files/ncov/open/reference/metadata.tsv.xz
      sequences: https://data.nextstrain.org/files/ncov/open/reference/sequences.fasta.xz
    - name: references
      metadata: data/references_metadata.tsv
      sequences: data/references_sequences.fasta

    builds:
      custom-build:
        title: "Build with custom data and example data"
        subsampling_scheme: all
        auspice_config: auspice-config-custom-data.json

    filter:
      ~{build_name}:
        min_length: 5000
        skip_diagnostics: True

    EOF
        export CONFIGFILE_FLAG="--configfile builds.yaml"
        echo "=====builds.yaml====="
        cat builds.yaml
        mv builds.yaml $INDIR/.
      fi

      cp ~{sequence_fasta} $INDIR/.
      cp ~{metadata_tsv} $INDIR/.
      wget https://raw.githubusercontent.com/nextstrain/ncov-tutorial/main/auspice-config-custom-data.json
      cat auspice-config-custom-data.json | sed 's/custom_data/~{build_name}/g' > ${INDIR}/auspice-config-custom-data.json

      echo "=====auspice-config-custom-data.json====="
      cat ${INDIR}/auspice-config-custom-data.json
    fi

    if [[ -n "~{configfile_yaml}" ]]; then
      export CONFIGFILE_FLAG="--configfile ~{configfile_yaml}"
    fi

    echo "CONFIGFILE_FLAG: " ${CONFIGFILE_FLAG}

    # If a tar gz of contextual sequences are provided such as GISAID Regional datasets, move it to the ncov folder
    if [[ -n "~{context_targz}" ]] ; then
      cp ~{context_targz} $INDIR/.
    fi

    # If a custom zipped folder of configs are provided, move it to ncov
    if [[ -n "~{custom_zip}" ]] ; then
      # Link custom profile (zipped version)
      cp ~{custom_zip} here_custom.zip
      CUSTOM_DIR=`unzip -Z1 here_custom.zip | head -n1 | sed 's:/::g'`
      unzip here_custom.zip
      cp -r $CUSTOM_DIR $INDIR/.

      # Draft: if passing config file from zip folder
      # BUILDYAML=`ls -1 $CUSTOM_DIR/*.yaml | head -n1`
      # cp $BUILDYAML $INDIR/build_custom.yaml
    fi

    # Max out the number of threads
    PROC=`nproc`

    # Run nextstrain
    nextstrain build \
      --cpus $PROC \
      --memory  ~{memory}Gib \
      --native $INDIR $CONFIGFILE_FLAG \
      ~{"--config active_builds=" + active_builds}

    if [[ -n "~{s3deploy}" ]] ; then
      # may be replaced with Nextstrain.org login instead, check docs
      # https://docs.nextstrain.org/projects/cli/en/latest/commands/remote/upload/#
      export AWS_ACCESS_KEY_ID=~{AWS_ACCESS_KEY_ID}
      export AWS_SECRET_ACCESS_KEY=~{AWS_SECRET_ACCESS_KEY}

      # deploy to Nextstrain Groups, todo: rename s3deploy
      nextstrain remote upload ~{s3deploy} $INDIR/auspice/*.json
    fi

    # Prepare output
    mv $INDIR/auspice .
    zip -r auspice.zip auspice

    # For debugging
    mv $INDIR/results .
    cp $INDIR/.snakemake/log/*.log results/.
    zip -r results.zip results
  >>>
  output {
    File auspice_zip = "auspice.zip"  # final output
    File results_zip = "results.zip"  # for debugging
  }
  runtime {
    docker: "nextstrain/base:latest"
    cpu : cpu
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
