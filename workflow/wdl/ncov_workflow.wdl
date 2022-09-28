version 1.0

import "tasks/nextstrain.wdl" as nextstrain

workflow Nextstrain_WRKFLW {
  input {
    # ncov
    # Option 1: Pass in a sequence and metadata files, create a configfile_yaml
    File? sequence_fasta
    File? metadata_tsv
    File? context_targz #<= optional contextual seqs in a tarball
    String? build_name
    
    # Option 2: Use a custom config file (e.g. builds.yaml) with https or s3 sequence or metadata files
    File? configfile_yaml
    File? custom_zip      # optional modifier: add a my_profiles.zip folder for my_auspice_config.json
    String? active_builds # optional modifier: specify "Wisconsin,Minnesota,Iowa"
    
    # Option 3? GISAID augur zip?
    # File? gisaid_zip # tarball
    
    # Optional Keys for deployment
    String? s3deploy
    String? AWS_ACCESS_KEY_ID
    String? AWS_SECRET_ACCESS_KEY
    
    # By default, run the ncov workflow (can swap it for zika or something else)
    String pathogen_giturl = "https://github.com/nextstrain/ncov/archive/refs/heads/master.zip"
    Int? cpu
    Int? memory       # in GiB
    Int? disk_size
  }

  call nextstrain.nextstrain_build as build {
    input:
      # Option 1
      sequence_fasta = sequence_fasta,
      metadata_tsv = metadata_tsv,
      context_targz = context_targz,
      build_name = build_name,
  
      # Option 2
      configfile_yaml = configfile_yaml,
      custom_zip = custom_zip,
      active_builds = active_builds,
  
      # Optional deploy to s3 site
      s3deploy = s3deploy,
      AWS_ACCESS_KEY_ID = AWS_ACCESS_KEY_ID,
      AWS_SECRET_ACCESS_KEY = AWS_SECRET_ACCESS_KEY,
  
      pathogen_giturl = pathogen_giturl,
      cpu = cpu,
      memory = memory,
      disk_size = disk_size
  }

  output {
    File auspice_zip = build.auspice_zip
    File results_zip = build.results_zip
  }
}
