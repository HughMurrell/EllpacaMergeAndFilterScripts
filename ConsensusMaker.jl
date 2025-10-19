using Pkg
Pkg.activate(".")
ENV["MPLBACKEND"] = "Agg"
using FASTX, BioSequences, StatsBase, DataFrames, CSV
# using SplitApplyCombine

########### first define some functions for later use

function degap(s::String)
    return replace(s,"-"=>"")
end

function write_fasta(filename, seqs;
    names=String[], LongSequence = false, append = false, aa = false)
    if !LongSequence
        if aa
            seqs = [BioSequences.LongAA(s) for s in seqs]
        else
            seqs = [BioSequences.LongDNA{4}(s) for s in seqs]
        end
    end
    stream = open(FASTA.Writer, filename, append=append)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, n) in zip(seqs, names)
        i += 1
        write(stream, FASTA.Record(n, s))
    end
    close(stream)
end

function read_fasta_with_names_and_descriptions(filename; seqtype=String)
    nams=[]
    seqs=[]
    reader = open(FASTA.Reader, filename)

        for rec in reader
            push!(nams,FASTA.description(rec))
            push!(seqs,FASTA.sequence(seqtype,rec))
        end

    close(reader)
    return nams, seqs
end

function read_fasta(filename; seqtype=String)
    nams=[]
    seqs=[]
    reader = open(FASTA.Reader, filename)

        for rec in reader
            push!(nams,FASTA.identifier(rec))
            push!(seqs,FASTA.sequence(seqtype,rec))
        end

    close(reader)
    return nams, seqs
end

function consensus(seqs)
    cons = join([mode([seqs[i][j]
                    for i in 1:length(seqs)])
                        for j in 1:length(seqs[1])])
    return(cons)
end

function agreement(ref,seq)
    if length(ref) == length(seq)
        return( sum(collect(ref).==collect(seq)) / length(ref) )
    else
        return 0.0
    end
end



############## end of function definitions

println("ConsensusMaker script")
println("using Julia version: $(VERSION)")
t1 = time()

nams, seqs = read_fasta_with_names_and_descriptions("panels/HIV1_ALL_2022_env_DNA.fasta")
@show length(seqs)
keeps = (x->startswith(x,"C.")).(nams)
@show sum(keeps)
write_fasta("panels/subtypeC_nt.fasta",seqs[keeps],names=nams[keeps])
cons=degap(consensus(seqs[keeps]))
write_fasta("panels/consensusC_env_nt.fasta",[cons],names=["consensusC"])
t2 = time()
println("All Renaming Done in $((t2 - t1)/60) minutes, see housekeeping_renaming.csv in reports/")
exit()

