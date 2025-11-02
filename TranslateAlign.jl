using Pkg
Pkg.activate(".")
ENV["MPLBACKEND"] = "Agg"
using FASTX, BioSequences, StatsBase, DataFrames, CSV
using MAFFT_jll
# using SplitApplyCombine

########### first define some functions for later use

function degap(s::String)
    return replace(s,"-"=>"")
end

function degap(s::LongDNA{4})
    return filter!(!isgap,s)
end

function degap(s::LongAA)
    return filter!(!isgap,s)
end

function my_mafft(seedpath, inpath, outpath)
    print("aligning $(basename(inpath))....")
    cmd = `$(mafft()) --quiet --thread 2 --ep 2 --op 3 --add $(inpath) --out $(outpath) $(seedpath) `
    run(cmd)
    println("completed!")
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

function read_fasta(filename; seqtype=LongDNA{4} )
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

function read_aa_fasta(filename; seqtype=LongAA )
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



############## end of function definitions ##############################

println("TranslateAlign script")
println("using Julia version: $(VERSION)")
t1 = time()
in_dir="ellpaca_collection/CAPRISA_002/env_nt_merged/"
out_dir="ellpaca_collection/CAPRISA_002/env_aa_aligned_hxb2/"
ref_file="panels/hxb2-env.fasta"
ref_file_aa="panels/hxb2-env-aa.fasta"
cons_file="panels/consensusC_env_nt.fasta"
cons_file_aa="panels/consensusC_env_aa.fasta"
seed_file_nt="panels/refcons_env_nt_ali.fasta"
seed_file_aa="panels/refcons_env_aa_ali.fasta"
mkpath(out_dir)
reject_dir="ellpaca_collection/CAPRISA_002/env_aa_aligned_hxb2_rejects/"
mkpath(reject_dir)

#println("loading reference sequence...")
#ref_nams, ref_seqs=read_fasta(ref_file)
#cons_nams, cons_seqs=read_fasta(cons_file)

# generate seed alignment
nnams, sseqs = read_fasta(ref_file)
sseqs=(BioSequences.translate).(sseqs)
write_fasta(ref_file_aa,sseqs,names=nnams,aa=true,append=false)
nnams, sseqs = read_fasta(cons_file)
sseqs=(BioSequences.translate).(sseqs)
write_fasta(cons_file_aa,sseqs,names=nnams,aa=true,append=false)
my_mafft(ref_file_aa,cons_file_aa,seed_file_aa)

println("translating sequences...")
file_names=readdir(in_dir)
file_names=file_names[(x->occursin(".fasta",x)).(file_names)]
hk=DataFrame( donor=[], sequences=[], filtered=[], no_start=[], bad_match=[] )
for file in file_names
    donor=split(file,"_")[2]
    in_file = in_dir * file
    out_file = out_dir * replace(file, "_nt_" => "_aa_")
    print("translating and " )
    nams, seqs = read_fasta(in_file)
    seq_count=length(seqs)
    # nams=vcat(ref_nams, cons_nams, nams)
    # seqs=vcat(ref_seqs, cons_seqs, (degap).(seqs))
    seqs=(BioSequences.translate).(seqs)
    
    write_fasta(out_file, seqs, names=nams, aa=true, append=false)
    ali_file=out_file[1:end-6] * "_a-maff_w-hxb2.fasta"
    my_mafft(seed_file_aa, out_file, ali_file)
    rm(out_file)
        
    # now trim ali_file to match ref
    ali_reject_file=reject_dir * file[1:end-6] * "_a-maff_w-hxb2+rejects.fasta"
    
    nams, seqs = read_aa_fasta(ali_file)
    # startTrim=min( findfirst('-'.!=(collect(string(seqs[1])))), findfirst('-'.!=(collect(string(seqs[2])))) )
    # stopTrim=max( findlast('-'.!=(collect(string(seqs[1])))), findlast('-'.!=(collect(string(seqs[2])))) )
    startTrim=findfirst('-'.!=(collect(string(seqs[1]))))
    stopTrim=findlast('-'.!=(collect(string(seqs[1]))))
    @show startTrim, stopTrim
    seqs=(x->x[startTrim:stopTrim]).(seqs)
    reject_seqs=[seqs[1],seqs[2]]
    reject_nams=[nams[1],nams[2]]
    
    # check for start
    keeps=(x->x[ findfirst('-'.!=(collect(string(x))))   ]==(AA_M)).(seqs)
    no_start_reject_count=sum((!).(keeps))
    reject_seqs=vcat(reject_seqs,seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*"_noStartReject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    seqs=seqs[keeps]
    
    # check for stop does not work
    # keeps=(x->x[end]!=(AA_Gap)).(seqs)
    # bad_term_reject_count=sum((!).(keeps))
    # reject_seqs=vcat(reject_seqs,seqs[(!).(keeps)])
    # reject_nams=vcat(reject_nams,(x->x*"_badTermReject").(nams[(!).(keeps)]))
    # nams=nams[keeps]
    # seqs=seqs[keeps]
    
    # check for minimum match with consensusC again
    match_ratios=(x->sum(collect(seqs[2]).==(collect(x)))/length(x)).(seqs)
    @show minimum(match_ratios)
    keeps=match_ratios.>0.8
    keeps[1]=true # dont reject the reference
    bad_match_count=sum((!).(keeps))
    reject_seqs=vcat(reject_seqs,seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*"_badMatch").(nams[(!).(keeps)]))
    nams=nams[keeps]
    seqs=seqs[keeps]
    
    write_fasta(ali_file,seqs[3:end], names=nams[3:end], aa=true, append=false)
    filtered_count=length(seqs)-2
    write_fasta(ali_reject_file,reject_seqs, names=reject_nams, aa=true, append=false)
    println("trimmed and filtered")
    
    #now realign
    my_mafft(seed_file_aa,ali_file,ali_file)
    println("and realigned")
    #house keeping
    hk_rec=[donor,seq_count,filtered_count,no_start_reject_count,bad_match_count]
    push!(hk,hk_rec)
end

# total the housekeeping records
rec = vcat(["totals"],sum.(eachcol(hk)[2:end]))
push!(hk, rec)

CSV.write(out_dir * "housekeeping_aa_maff_filter.csv",hk)

t2 = time()
println()
println("All files translated and aligned in $((t2 - t1)/60) minutes.")

