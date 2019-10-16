"""
This julia script takes an input CCS.fastq file along with a primer sheet specifying
the PCR design and demultiplexes into individual .fastq files. For usage with
Snakemake pipeline.
"""

input_fastq_path = snakemake.input[1]
dataset = split(basename(input_fastq_path),".")[1]
error_rate = 0.01
templates = snakemake.params["templates"]

println("Running demux with the following parameters...")
println("Input fastq: $(input_fastq_path)")
#println("Templates")
#for template in templates
#    println(template)
#end
println("Target size: $(snakemake.params["target_size"])")
println("Index type: $(snakemake.params["index_type"])")
println("Error rate: $(error_rate)")

#load required packages
ENV["MPLBACKEND"] = "Agg"
using NextGenSeqUtils, DataFrames, DataFramesMeta, CSV, StatsBase, IterTools

#defining functions
function unique_not_substr(a)
    out = []
    for i in unique(a)
        res = true
        for j in unique(a)
            if occursin(i, j) & (i != j)
                res = false
            end
        end
        if res
            push!(out, i)
        end
    end
    return out
end

function iterative_primer_match(seqs,full_primers,window::Int,slide_by::Int;tol_one_error=true)
    if(slide_by + window - 1 > minimum(length.(full_primers)))
        @warn ("Matching window extends beyond shortest primer. This is ok, but check that you aren't matching something too short.")
    end
    primers = [p[1:min(window,minimum(length.(full_primers)))] for p in full_primers]
    filter = fast_primer_match(seqs,primers,tol_one_error=tol_one_error);
    for i in 2:slide_by
        unresolved = filter .== 0
        primers = [p[i:min(i + window - 1,minimum(length.(full_primers)))] for p in full_primers]
        filter[unresolved] = fast_primer_match(seqs[unresolved],primers,tol_one_error=tol_one_error);
    end
    return filter
end

function sliding_demux_dict(seqs,fwd_primers,window::Int,slide_by::Int; verbose = true, phreds = nothing, tol_one_error = true)
    fwd_matches = iterative_primer_match(seqs,fwd_primers,window,slide_by,tol_one_error=tol_one_error)
    rev_comp_bool = fwd_matches .< 0
    keepers = abs.(fwd_matches) .> 0
    fwd_matches = abs.(fwd_matches)
    pair_keeps = fwd_matches[keepers]
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            println(s[1], " => ", s[2])
        end
    end
    if phreds == nothing
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]

                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            end
        end
        return seq_dict
    else
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Vector{Int8},Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            end
        end
        return seq_dict
    end
end

#define nextera indexes
N7_dic = Dict(
    "N701" => "TCGCCTTA",
    "N702" => "CTAGTACG",
    "N703" => "TTCTGCCT",
    "N704" => "GCTCAGGA",
    "N705" => "AGGAGTCC",
    "N706" => "CATGCCTA",
    "N707" => "GTAGAGAG",
    "N708" => "CCTCTCTG",
    "N709" => "AGCGTAGC",
    "N710" => "CAGCCTCG",
    "N711" => "TGCCTCTT",
    "N712" => "TCCTCTAC"
);

S5_dic = Dict(
    "S501" => "TAGATCGC",
    "S502" => "CTCTCTAT",
    "S503" => "TATCCTCT",
    "S504" => "AGAGTAGA",
    "S505" => "GTAAGGAG",
    "S506" => "ACTGCATA",
    "S507" => "AAGGAGTA",
    "S508" => "CTAAGCCT",
    "S517" => "GCGTAAGA"
);

#define universal adapter sequences
N7_univ = "CAAGCAGAAGACGGCATACGAGAT";
S5_univ = "AATGATACGGCGACCACCGAGATCTACAC";

N7_suffix = "GTCTCGTGGGCTCGG"
S5_suffix = "TCGTCGGCAGCGTC"

"""
demux CCS based on nextera illumina adapter sequences using a sliding primer match.
Writes collections of reads to .fastq named by index.
"""
function demux_nextera(file; outdir = "demux/", thresh = 100, verbose = true)
    if verbose println("Demultiplexing $(file)...") end
    outdir = strip(outdir,'/')
    seqs, phreds, names = read_fastq(file);
    if verbose println("Finding sequences with N7 index...") end
    #sliding window demultiplex on forward primers
    fwd_demux_dic = sliding_demux_dict(seqs,
                                       [N7_univ],
                                       10,
                                       14,
                                       verbose=false,
                                       phreds = phreds)
    #There is only one entry in this dict
    #retrieve from demux_dic
    seqs_N7 = [i[1] for i in fwd_demux_dic[1]];
    phreds_N7 = [i[2] for i in fwd_demux_dic[1]];
    names_N7 = names[[i[3] for i in fwd_demux_dic[1]]];
    if verbose println("$(length(seqs_N7)) sequences match N7 adapters.") end
    #match to reverse adapter
    S5_matches = iterative_primer_match(seqs_N7, [S5_univ], 10, 19, tol_one_error=true);
    S5_keepers = S5_matches .< 0;
    #filter to reverse adapter matches
    seqs_N7S5 = seqs_N7[S5_keepers];
    phreds_N7S5 = phreds_N7[S5_keepers];
    names_N7S5 = names_N7[S5_keepers];
    if verbose println("$(length(seqs_N7S5)) sequences match N7 and S5 adapters.") end
    seqs_N7S5_trim = [double_primer_trim(s, p, N7_univ, S5_univ) for (s,p) in zip(seqs_N7S5, phreds_N7S5)];
    #run demux by paired indexes
    nextera_demux_dic = demux_dict(
        [s for (s,p) in seqs_N7S5_trim],
        collect(values(N7_dic)),
        collect(values(S5_dic));
        phreds = [p for (s,p) in seqs_N7S5_trim],
        tol_one_error = false,
        verbose = false);
    return nextera_demux_dic, names_N7S5
end

t1 = time()
#filter .fastq
#filtered_path = snakemake.output[2] #julia temp paths don't work
#println("Filtering .fastq file...")
#@time fastq_filter(input_fastq_path,
#                   filtered_path, #path here
#                   error_rate = error_rate,
#                   min_length = snakemake.params["target_size"]*0.8,
#                   max_length = snakemake.params["target_size"]*1.2)

if snakemake.params["index_type"] == "nextera"
    nextera_demux_dic,seqnames = demux_nextera(input_fastq_path, outdir = "demux/$(dataset)/nextera/")
    nex_tuples = collect(keys(nextera_demux_dic))
    N7 = collect(keys(N7_dic))
    S5 = collect(keys(S5_dic))
    index_tuples = [(N7[x],S5[y]) for (x,y) in nex_tuples]
    indexes2tuples = Dict(zip(index_tuples,nex_tuples))
    for template in collect(keys(templates))
        indexes = templates[template]
        if (indexes["Index_1"],indexes["Index_2"]) in index_tuples
            template_seqs = nextera_demux_dic[indexes2tuples[(indexes["Index_1"],indexes["Index_2"])]]
            if length(template_seqs) < 10 @warn "Less than 10 reads for $(template): $(indexes)" end
            write_fastq(snakemake.output[1]*"/$(template).fastq",
                        [i[1] for i in template_seqs],
                        [i[2] for i in template_seqs];
                        names = seqnames[[i[3] for i in template_seqs]])
        else
            @warn "No reads found for $(template): $(indexes)"
        end
    end
else
    @warn "index_type $(snakemake.params["index_type"]) not recognized."
end
t2 = time()
println("Demultiplex took $(t2 - t1) seconds.")
