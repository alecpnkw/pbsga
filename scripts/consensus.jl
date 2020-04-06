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
cons_collection, name_collection = denoiseDir(ARGS[1])
write_fasta(ARGS[2],
    cons_collection;
    names = name_collection)
t2 = time()
println("Consensus generation took $(t2-t1) seconds.")
