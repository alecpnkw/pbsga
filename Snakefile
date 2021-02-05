configfile: "config.yaml"

rule all:
    input:
        expand("consensus/{dataset}.fasta", dataset = config["datasets"])

rule demux:
    input:
        "raw-reads/{dataset}.fastq"
    output:
        directory("demux/{dataset}/"),
        counts = "counts/{dataset}.csv",
        filtered_file = temporary("tmp/{dataset}_filt.fastq"),
    params:
        index_type = lambda wc: config["datasets"][wc.dataset]["index_type"],
        outdir = "demux/{dataset}",
        templates = lambda wc: config["datasets"][wc.dataset]["templates"]
    script:
        "scripts/demux.jl"

rule consensus:
    input:
        "demux/{dataset}/"
    output:
        "consensus/{dataset}.fasta"
    threads: 6
    script:
        "scripts/consensus.jl"
