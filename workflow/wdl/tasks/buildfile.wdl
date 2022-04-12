version 1.0

task mk_buildconfig {
  input {
    File sequence_fasta  # Could change this to Array[File] and loop the HEREDOC
    File metadata_tsv
    String build = "example"
    String dockerImage = "nextstrain/base:latest"
    Int cpu = 1
    Int disk_size = 5
    Float memory = 3.5
  }
  command <<<
    cat << EOF > build.yaml
    inputs:
    - name: ~{build}
      metadata: ~{basename(metadata_tsv)}
      sequences: ~{basename(sequence_fasta)}
    - name: references
      metadata: data/references_metadata.tsv
      sequences: data/references_sequences.fasta
    EOF
  >>>
  output {
    File buildconfig = "build.yaml"
  }
  runtime {
    docker: dockerImage
    cpu : cpu
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

# Public Reference Datasets in case we want to add default "context" strains for ncov
# nextstrain remote ls s3://nextstrain-data/files/ncov/open &> list.txt
# files/ncov/open/LOCATION/metadata.tsv.xz where LOCATION = one of ['africa', 'asia', 'europe', 'global', 'north-america', 'oceania', 'south-america']
# files/ncov/open/LOCATION/sequences.fasta.xz
# TODO: since these are all s3, just generate string (avaid download and localazation/delocalization) for the mk_buildconfig
# Clarification: Do these need to be sampled down?