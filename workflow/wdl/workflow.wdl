version 1.0

# import "tasks/nextstrain.wdl" as nextstrain # <= modular method
import "tasks/buildfile.wdl" as buildfile
import "tasks/nextstrain.wdl" as nextstrain
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

  if (defined(sequence_fasta)) {
    call buildfile.mk_buildconfig as mk_buildconfig {
      input:
        sequence_fasta = select_first([sequence_fasta]),
        metadata_tsv = select_first([metadata_tsv]),
        build = build_name,
        dockerImage = docker_path
    }
  }

  # call nextstrain.nextstrain build as build {  # <= modular method
  call nextstrain.nextstrain_build as build {
    input:
      sequence_fasta = sequence_fasta,
      metadata_tsv = metadata_tsv,
      build_yaml = select_first([build_yaml, mk_buildconfig.buildconfig]), # Accepts Option 1 or Option 2
      custom_zip = custom_zip,
      cpu = cpu,
      memory = memory,
      disk_size = disk_size,
      dockerImage = docker_path,
      pathogen_giturl = pathogen_giturl,
      active_builds = active_builds,
      AWS_ACCESS_KEY_ID = AWS_ACCESS_KEY_ID,
      AWS_SECRET_ACCESS_KEY = AWS_SECRET_ACCESS_KEY,
      s3deploy = s3deploy
  }

  output {
    #Array[File] json_files = build.json_files
    File auspice_zip = build.auspice_zip
    File results_zip = build.results_zip
  }
}
