version 1.0

# import "tasks/buildfile.wdl" as buildfile
import "tasks/nextstrain.wdl" as nextstrain  # <= modular method
# import "tasks/ncov_ingest.wdl" as ncov_ingest

workflow Nextstrain_WRKFLW {
  input {
    # Option 1: Pass in a sequence and metadata files, create a build_yaml
    File? sequence_fasta
    File? metadata_tsv
    String? build_name

    # Option 2: Use a custom build_yaml file with https or s3 sequence or metadata files
    File? build_yaml
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
    String docker_path = "nextstrain/base:latest"
    Int? cpu
    Int? memory       # in GiB
    Int? disk_size
  }

  call nextstrain.nextstrain_build as build {
    input:
      # Option 1
      sequence_fasta = sequence_fasta,
      metadata_tsv = metadata_tsv,
      build_name = build_name,

      # Option 2
      build_yaml = build_yaml,
      custom_zip = custom_zip,
      active_builds = active_builds,

      # Optional deploy to s3 site
      s3deploy = s3deploy,
      AWS_ACCESS_KEY_ID = AWS_ACCESS_KEY_ID,
      AWS_SECRET_ACCESS_KEY = AWS_SECRET_ACCESS_KEY,

      pathogen_giturl = pathogen_giturl,
      dockerImage = docker_path,
      cpu = cpu,
      memory = memory,
      disk_size = disk_size
  }

  output {
    #Array[File] json_files = build.json_files
    File auspice_zip = build.auspice_zip
    File results_zip = build.results_zip
  }
}
