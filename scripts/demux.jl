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
#println("Target size: $(snakemake.params["target_size"])")
println("Index type: $(snakemake.params["index_type"])")
println("Error rate: $(error_rate)")

#load required packages
ENV["MPLBACKEND"] = "Agg"
using NextGenSeqUtils, DataFrames, DataFramesMeta, CSV,
StatsBase, IterTools, StringDistances, Statistics

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

#"sga" primers
sga_primers_f = Dict(
  "SGA_F01" => "CTACACTCGCCTTATCGTCGGCAGCGTC",
  "SGA_F02" => "CTACACCTAGTACGTCGTCGGCAGCGTC",
  "SGA_F03" => "CTACACTTCTGCCTTCGTCGGCAGCGTC",
  "SGA_F04" => "CTACACGCTCAGGATCGTCGGCAGCGTC",
  "SGA_F05" => "CTACACAGGAGTCCTCGTCGGCAGCGTC",
  "SGA_F06" => "CTACACCATGCCTATCGTCGGCAGCGTC",
  "SGA_F07" => "CTACACGTAGAGAGTCGTCGGCAGCGTC",
  "SGA_F08" => "CTACACCAGCCTCGTCGTCGGCAGCGTC",
  "SGA_F09" => "CTACACTGCCTCTTTCGTCGGCAGCGTC",
  "SGA_F10" => "CTACACTCCTCTACTCGTCGGCAGCGTC",
  "SGA_F11" => "CTACACTCATGAGCTCGTCGGCAGCGTC",
  "SGA_F12" => "CTACACCCTGAGATTCGTCGGCAGCGTC",
  "SGA_F13" => "CTACACTAGCGAGTTCGTCGGCAGCGTC",
  "SGA_F14" => "CTACACGTAGCTCCTCGTCGGCAGCGTC",
  "SGA_F15" => "CTACACTACTACGCTCGTCGGCAGCGTC",
  "SGA_F16" => "CTACACAGGCTCCGTCGTCGGCAGCGTC",
  "SGA_F17" => "CTACACGCAGCGTATCGTCGGCAGCGTC",
  "SGA_F18" => "CTACACCTGCGCATTCGTCGGCAGCGTC",
  "SGA_F19" => "CTACACGAGCGCTATCGTCGGCAGCGTC",
  "SGA_F20" => "CTACACCGCTCAGTTCGTCGGCAGCGTC",
  "SGA_F21" => "CTACACGTCTTAGGTCGTCGGCAGCGTC",
  "SGA_F22" => "CTACACACTGATCGTCGTCGGCAGCGTC",
  "SGA_F23" => "CTACACTAGCTGCATCGTCGGCAGCGTC",
  "SGA_F24" => "CTACACGACGTCGATCGTCGGCAGCGTC"
);

sga_primers_r = Dict(
  "SGA_R01" => "CGAGATCTCTCTATGTCTCGTGGGCTCGG",
  "SGA_R02" => "CGAGATTATCCTCTGTCTCGTGGGCTCGG",
  "SGA_R03" => "CGAGATGTAAGGAGGTCTCGTGGGCTCGG",
  "SGA_R04" => "CGAGATACTGCATAGTCTCGTGGGCTCGG",
);

SGA_F_univ = "CTACACNNNNNNNNTCGTCGGCAGCGTC"
SGA_R_univ = "CGAGATNNNNNNNNGTCTCGTGGGCTCGG"

function find_nextera_suffix(query_seq, query_phred, suffix; start_ix = 34, end_ix = 12, try_reverse_comp = true)
    for ix in start_ix:-1:end_ix
        window = query_seq[ix : ix + length(suffix) - 1]
        d = evaluate(Hamming(), window, suffix)
        if d < 2
            return(query_seq[ix - 8:end], query_phred[ix - 8:end])
        end
    end
    #try reverse complement
    if try_reverse_comp
        query_seq = reverse_complement(query_seq)
        query_phred = query_phred[end:-1:1]
        for ix in start_ix:-1:end_ix
            window = query_seq[ix : ix + length(suffix) - 1]
            d = evaluate(Hamming(), window, suffix)
            if d < 2
                return(query_seq[ix - 8:end], query_phred[ix - 8:end])
            end
        end
    end
end

function get_nextera_matches(seqs, phreds)
    #define universal adapter sequences
    N7_univ = "CAAGCAGAAGACGGCATACGAGAT";
    S5_univ = "AATGATACGGCGACCACCGAGATCTACAC";

    N7_suffix = "GTCTCGTGGGCTCGG";
    S5_suffix = "TCGTCGGCAGCGTC";

    #find N7
    N7_matches = find_nextera_suffix.(seqs, phreds, N7_suffix;
        start_ix = length(N7_univ) + 10, end_ix = 12, try_reverse_comp = true)

    N7_coords = [i for (i,m) in enumerate(N7_matches) if !isnothing(m)]
    N7_keeps = N7_matches[N7_coords]

    #find S5 (no revc)
    matches = find_nextera_suffix.([s for (s,p) in N7_keeps],
        [p for (s,p) in N7_keeps], S5_suffix;
        start_ix = length(S5_univ) + 10, end_ix = 12, try_reverse_comp = true) #find a better solution than this
    coords = [i for (i,m) in enumerate(matches) if !isnothing(m)]
    keeps = matches[coords]

    return [s for (s,p) in keeps], [p for (s,p) in keeps], coords
end

"""
demux CCS based on nextera illumina adapter sequences using a sliding primer match.
Writes collections of reads to .fastq named by index.
"""
function demux_nextera(file; verbose = true)
    if verbose println("Demultiplexing $(file)...") end
    seqs, phreds, seqnames = read_fastq(file);

    #proper usage
    @time matched_seqs, matched_phreds, coords = get_nextera_matches(seqs, phreds);
    names_N7S5 = seqnames[coords];
    nextera_demux_dic = demux_dict(matched_seqs,collect(values(N7_dic)),collect(values(S5_dic));
        phreds = matched_phreds,tol_one_error = false,verbose = false);
    return nextera_demux_dic, names_N7S5
end

t1 = time()
#filter .fastq
filtered_path = snakemake.output[2] #julia temp paths don't work
println("Filtering .fastq file...")
@time fastq_filter(input_fastq_path,
                   filtered_path, #path here
                   error_rate = error_rate,
                   min_length = 400,
                   max_length = 5000)

if snakemake.params["index_type"] == "nextera"
    nextera_demux_dic,seqnames = demux_nextera(filtered_path)
    nex_tuples = collect(keys(nextera_demux_dic))
    N7 = collect(keys(N7_dic))
    S5 = collect(keys(S5_dic))
    index_tuples = [(N7[x],S5[y]) for (x,y) in nex_tuples]
    indexes2tuples = Dict(zip(index_tuples,nex_tuples))
    for template in collect(keys(templates))
        indexes = templates[template]
        if (indexes["N7_Index"],indexes["S5_Index"]) in index_tuples
            template_seqs = nextera_demux_dic[indexes2tuples[(indexes["N7_Index"],indexes["S5_Index"])]]
            #match template sequences, length here

            if length(template_seqs) < 3 @warn "Less than 3 reads for $(template): $(indexes)" end
            trimmed_seqs = [
                double_primer_trim(s,p,
                N7_suffix*templates[template]["Reverse_Primer_2ndRd_Sequence"],S5_suffix*templates[template]["Forward_Primer_2ndRd_Sequence"];
                buffer = 8)
            for (s,p) in template_seqs
            ]
            write_fastq(snakemake.output[1]*"/$(template).fastq",
                        [i[1] for i in trimmed_seqs],
                        [i[2] for i in trimmed_seqs];
                        names = seqnames[[i[3] for i in template_seqs]])
        else
            @warn "No reads found for $(template): $(indexes)"
        end
    end
elseif snakemake.params["index_type"] == "sga_primer"
    seqs, phreds, seqnames = read_fastq(filtered_path)
    demux_dic = demux_dict(seqs,
        [i[1:16] for i in collect(values(sga_primers_f))],
        [i[1:16] for i in collect(values(sga_primers_r))];
        phreds = phreds,
        tol_one_error = true);
    nex_tuples = collect(keys(demux_dic))
    SGA_F = collect(keys(sga_primers_f))
    SGA_R = collect(keys(sga_primers_r))
    index_tuples = [(SGA_F[x],SGA_R[y]) for (x,y) in nex_tuples]
    indexes2tuples = Dict(zip(index_tuples,nex_tuples))
    for template in collect(keys(templates))
        indexes = templates[template]
        if (indexes["N7_Index"],indexes["S5_Index"]) in index_tuples
            template_seqs = demux_dic[indexes2tuples[(indexes["N7_Index"],indexes["S5_Index"])]]
            println(length(template_seqs))
            index_trimmed = [double_primer_trim(s,p,SGA_F_univ,SGA_R_univ) for (s,p,n) in template_seqs];
            println(length(index_trimmed))
            #match template
            keeps = iterative_primer_match([s for (s,p) in index_trimmed], [templates[template]["Forward_Primer_2ndRd_Sequence"]],12,25; tol_one_error=true) .> 0 #primer matching, primers needs to be array
            seqs_keeping = index_trimmed[keeps]
            if length(seqs_keeping) == 0 @warn "No reads for $(template): $(indexes)"; continue end
            println(length(seqs_keeping))
            #length filter
            filtered_seqs = length_filter(
                [s for (s,p) in seqs_keeping],
                [p for (s,p) in seqs_keeping],
                seqnames[[i[3] for i in template_seqs]][keeps],
                Int(round(median(length.([s for (s,p) in seqs_keeping]))*0.9)),
                Int(round(median(length.([s for (s,p) in seqs_keeping]))*1.1))
                ) #check
            filtered_seqs = collect(zip(filtered_seqs...))

            if length(filtered_seqs) < 3
                @warn "Less than 3 reads for $(template): $(indexes)"
                continue
            end
            trimmed_seqs = [
                double_primer_trim(s,p,
                templates[template]["Forward_Primer_2ndRd_Sequence"],templates[template]["Reverse_Primer_2ndRd_Sequence"])
            for (s,p,n) in filtered_seqs
            ]
            write_fastq(snakemake.output[1]*"/$(template).fastq",
                        [i[1] for i in trimmed_seqs],
                        [i[2] for i in trimmed_seqs];
                        names = [n for (s,p,n) in filtered_seqs])
        else
            @warn "No reads found for $(template): $(indexes)"
        end
    end
else
    @warn "index_type $(snakemake.params["index_type"]) not recognized."
end
t2 = time()
println("Demultiplex took $(t2 - t1) seconds.")
