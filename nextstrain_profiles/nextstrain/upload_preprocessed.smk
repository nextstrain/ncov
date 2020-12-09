localrules: download_metadata, download_sequence, download, upload, clean

rule upload:
    message: "Uploading intermediate files to {params.s3_bucket}"
    input:
        "results/masked.fasta",
        "results/aligned.fasta",
        "results/filtered.fasta",
        "results/sequence-diagnostics.tsv",
        "results/flagged-sequences.tsv",
        "results/to-exclude.txt"
    params:
        s3_bucket = config["S3_BUCKET"],
        compression = config["preprocess"]["compression"]
    run:
        for fname in input:
            shell(f"./scripts/upload-to-s3 {fname} s3://{params.s3_bucket}/{os.path.basename(fname)}.{params.compression}")
