using Distributed
@everywhere ENV["MPLBACKEND"] = "Agg"
@everywhere using RobustAmpliconDenoising, NextGenSeqUtils

function generateConsensusFromDir(dir, template_name)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        exit()
    end
    cons_collection = pmap(ConsensusFromFastq, files)
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    return seq_collection, seqname_collection
end

@everywhere function ConsensusFromFastq(file)
    seqs,phreds,seq_names = read_fastq(file)
    draft = consensus_seq(seqs)
    draft2 = refine_ref(draft, seqs)
    final_cons = refine_ref(draft2,seqs)
    alignments, maps, matches, matchContent = getReadMatches(final_cons, seqs, 0)
    cons_name = split(basename(file),".")[1]*" num_CCS=$(length(seqs)) min_agreement=$(round(minimum(matches); digits = 2))"
    return final_cons, cons_name
end

@everywhere begin
"""
Returns an array of degapped coordinates, such that
coords(ref, read)[i] gives you the position the aligned read/ref
that matches the i'th ungapped position in ref.
"""
function coords(ref, read)
    if length(ref) != length(read)
        error("Aligned strings are meant to be the same length.")
    end
    degappedRef = degap(ref)
    coordMap = zeros(Int64, length(degappedRef))
    count = 1
    for i in 1:length(degappedRef)
        while ref[count] == '-'
            count += 1
        end
        coordMap[i] = count
        count += 1
    end
    return coordMap
end
end

@everywhere begin
"""
Return matches to a candidate reference from a set of reads.
"""
function getReadMatches(candidate_ref, reads, shift; degap_param = true, kmer_align = true)
    alignments = []
    if kmer_align
        alignments = map(i -> kmer_seeded_align(candidate_ref, i), reads)
    else
        alignments = map(i -> nw_align(candidate_ref, i), reads)
    end

    maps = [coords(i...) for i in alignments]

    if (degap_param)
        matchContent = [[degap(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
        matches = [freq(matchContent[k], degap(candidate_ref[k:k+shift])) for k in 1:length(matchContent)]
    else
       matchContent = [[(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
       matches = [freq(matchContent[k], candidate_ref[k:k+shift]) for k in 1:length(matchContent)]
    end
    return alignments, maps, matches, matchContent
end
end

#load required packages
using Distributed
@everywhere ENV["MPLBACKEND"] = "Agg"
@everywhere using NextGenSeqUtils, RobustAmpliconDenoising

function denoiseDir(dir; fine_radius = 1.0)
    dir = strip(dir,'/')
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) .fastq files")
    else
        println("WARNING: no input .fastq files")
        exit()
    end
    cluster_collection = pmap(x -> threshDenoiseFASTQ(x; fine_radius = fine_radius), files)
    cons_collection = vcat([s for (s,n) in cluster_collection]...)
    name_collection = vcat([n for (s,n) in cluster_collection]...)
    return cons_collection, name_collection
end

@everywhere function threshDenoiseFASTQ(file; fine_radius = 1.0, thresh = 0.65)
    seqs,phreds,seq_names = read_fastq(file)
    if length(seqs) < 10
        println("Insufficient read depth for $(basename(file))!")
        return ([],[])
    end
    cluster_seqs, sizes, members = denoise(seqs; fine_radius = fine_radius)
    if length(cluster_seqs) < 1
        println("No viable clusters for $(basename(file))!")
        return ([],[])
    else
    props = round.(sizes ./ sum(sizes), digits = 2)
    cluster_names = [split(basename(file),'.')[1]*"_$(sizes[i])_$(props[i])" for i in 1:length(sizes)]
    end
    return cluster_seqs, cluster_names
end

#Calculate consensus sequences for each family.
t1 = time()
cons_collection, name_collection = generateConsensusFromDir(ARGS[1], "")
write_fasta(ARGS[2],
    cons_collection;
    names = name_collection)
t2 = time()
println("Consensus generation took $(t2-t1) seconds.")
