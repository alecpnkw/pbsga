from snakemake.remote.SFTP import RemoteProvider
SFTP = RemoteProvider(username="appankow", private_key="/home/appankow/.ssh/id_rsa_themis")

configfile: "AMP-04-config.yaml"

rule all:
    input:
        expand("consensus/{dataset}.fasta", dataset = config["datasets"])

rule demux:
    input:
        SFTP.remote("hercules/opt/shared/PacBio_PipelineData/{dataset}/{dataset}.fastq")
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
    shell:
        "julia -p {threads} scripts/consensus.jl {input} {output}"
