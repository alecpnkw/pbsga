configfile: "config.yaml"

rule all:
    input:
        expand("consensus/{dataset}.fasta", dataset = config["datasets"])

rule demux:
    input:
        "raw_reads/{dataset}.fastq"
    output:
        directory("demux/{dataset}/"),
        #temporary("tmp/{dataset}_filt.fastq")
    params:
        index_type = lambda wc: config["datasets"][wc.dataset]["index_type"],
        target_size = lambda wc: config["datasets"][wc.dataset]["target_size"],
        templates = lambda wc: config["datasets"][wc.dataset]["templates"]
    script:
        "scripts/demux.jl"

rule consensus:
    input:
        "demux/{dataset}/"
    output:
        "consensus/{dataset}.fasta"
    threads: 6
    shell:
        "julia -p {threads} scripts/consensus.jl {input} {output}"
