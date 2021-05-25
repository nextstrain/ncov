# See docs/multiple_inputs.md for an explanation of this file!

rule make_starting_files:
    message:
        """
        Creating starting files for the multiple inputs tutorial by decompressing {input.archive}
        """
    input:
        archive = "data/example_multiple_inputs.tar.xz"
    output:
        # Note: the command doesn't use these, but adding them here makes snakemake
        # aware that this rule produces them
        aus_meta = "data/example_metadata_aus.tsv",
        aus_seqs = "data/example_sequences_aus.fasta",
        world_meta = "data/example_metadata_worldwide.tsv",
        world_seqs = "data/example_sequences_worldwide.fasta"
    shell:
        """
        tar xf {input.archive} --directory data/
        """
