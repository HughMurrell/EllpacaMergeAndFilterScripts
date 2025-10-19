using Pkg
Pkg.activate(".")
ENV["MPLBACKEND"] = "Agg"
using FASTX, BioSequences, StatsBase, DataFrames, CSV
using SeededAlignment

########### first define some functions for later use

function degap(s::String)
    return replace(s,"-"=>"")
end

function degap(s::LongDNA{4})
    return filter!(!isgap,s)
end
        

function translate_to_aa(s::String)
    s=s[1:3*div(length(s),3)]
    rna = convert(LongRNA{4}, LongDNA{4}(s))
    return string(translate(rna))
end

function reverse_translate_with_gaps(aa_seq::String, nu_seq::String)
    nu_seq = degap(nu_seq)
    if 3*length(degap(aa_seq)) != length(nu_seq)
        # println("3 * length of degapped aa_seq = $(3*length(degap(aa_seq)))")
        # println(" length of nu_seq = $(length(nu_seq) )")
        # println("reading frame error")
        return("")
    end
    nu_seq = reverse(collect(nu_seq))
    ret = []
    for ch in collect(aa_seq)
        if ch == '-'
            push!(ret,'-')
            push!(ret,'-')
            push!(ret,'-')
        else
            push!(ret,pop!(nu_seq))
            push!(ret,pop!(nu_seq))
            push!(ret,pop!(nu_seq))
        end
    end
    return(join(ret))
end
        


"""
    my_mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
Julia wrapper for mafft.
"""
function my_mafft(inpath, outpath)
    # cmd = `mafft-fftns --quiet --thread 2 --ep 2 --op 3 --out $outpath $inpath`
    cmd = `mafft --quiet --thread 2 --ep 2 --op 3 --out $outpath $inpath`
    println(cmd)
    run(cmd)
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


function all_in_one_filter(in_dir,in_file,ref_nam,ref_seq,out_dir;hk=DataFrame())
    nams, seqs = read_fasta(in_dir * in_file)
    @show in_file, length(seqs)
    nams=vcat([ref_nam],nams)
    seqs=vcat([ref_seq],seqs)
    ali_seqs = msa_codon_align(ref_seq, seqs)
    startTrim=findfirst('-'.!=(collect(string(ali_seqs[1]))))
    stopTrim=findlast('-'.!=(collect(string(ali_seqs[1]))))
    @show startTrim, stopTrim
    trim_ali_seqs=(x->x[startTrim:stopTrim]).(ali_seqs)
    write_fasta(out_dir*in_file[1:end-6]*"_f-hmff.fasta",trim_ali_seqs,seq_names=nams)
end

function longest_open_reading_frame(cons)
    cons=degap(cons)
    cons_frames = [ cons[1:3*div(length(cons),3)],
                    cons[2:3*div(length(cons),3)-2],
                    cons[3:3*div(length(cons),3)-1] ]
    # then translate them
    cons_aa_frames = (x->String(BioSequences.translate(LongDNA{4}(degap(x))))).(cons_frames)
    cons_aa = cons_aa_frames[1]
        
    # find longest coding region in all frames
    best_reading_frame = 1
    start=1
    stop=length(cons_aa)
    max_match_length=0
    for rf in 1:3
        cons_aa = cons_aa_frames[rf]
        m=match(r"M[^\*]*\*",cons_aa)
        while ! isnothing(m)
            if length(m.match)>max_match_length
                best_reading_frame=rf
                start=m.offset
                stop=start-1+length(m.match)
                max_match_length=length(m.match)
            end
            m=match(r"M[^\*]*\*",cons_aa,m.offset+1)
        end
    end
    cons_aa = cons_aa_frames[best_reading_frame]
    cons_aa_trim = cons_aa[start:stop]
    cons_trim = cons[((start-1)*3+best_reading_frame):((stop)*3+best_reading_frame-1)]
    return cons_trim, cons_aa_trim, start, stop
end


function pairwise_filter(in_dir,in_file,ref_nam,ref_seq,out_dir,out_reject_dir;hk=DataFrame())
    score_params = ScoringScheme(edge_ext_begin= true, edge_ext_end = true )
    # seed_chain_align(A,B scoring=score_params)
    donor=split(in_file,"_")[2]
    visit=split(in_file,"_")[3]
    pool=split(split(in_file,"_")[6],"-")[3][2:end]
    facility=split(split(in_file,"_")[6],"-")[4][1:end]
    nams, seqs = read_fasta(in_dir * in_file)
    start_count = length(seqs)
    for i in 1:length(seqs)
        seq_trim, seq_aa_trim, start, stop = longest_open_reading_frame(seqs[i])
        pw_align = seed_chain_align(ref=ref_seq, query=seq_trim, scoring=score_params)
        startTrim=findfirst('-'.!=(collect(string(pw_align[1]))))
        stopTrim=findlast('-'.!=(collect(string(pw_align[1]))))
        seqs[i]=filter!(!isgap, pw_align[2][startTrim:stopTrim])
    end
    reject_seqs=[]
    reject_nams=[]
    keeps=(x->count(==(DNA_N),collect(x))==0).(seqs)
    orf_reject_count=start_count-sum(keeps)
    reject_seqs=vcat(reject_seqs,seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*"_orfReject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    seqs=seqs[keeps]
    nams=vcat([ref_nam],nams)
    ali_seqs = msa_codon_align(ref_seq, seqs, scoring=score_params)
    startTrim=findfirst('-'.!=(collect(string(ali_seqs[1]))))
    stopTrim=findlast('-'.!=(collect(string(ali_seqs[1]))))
    trim_ali_seqs=(x->x[startTrim:stopTrim]).(ali_seqs)
    keeps=(x->(x[1:3]==dna"ATG")).(trim_ali_seqs)
    no_start_codon_count=length(trim_ali_seqs)-sum(keeps)
    reject_seqs=vcat(reject_seqs,trim_ali_seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*"_noStartReject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    trim_ali_seqs=trim_ali_seqs[keeps]
    keeps=(x->(x[end-2:end]!=dna"---")).(trim_ali_seqs)
    no_stop_codon_count=length(trim_ali_seqs)-sum(keeps)
    reject_seqs=vcat(reject_seqs,trim_ali_seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*"_noStopReject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    trim_ali_seqs=trim_ali_seqs[keeps]
    match_ratios=(x->sum(collect(trim_ali_seqs[1]).==(collect(x)))/length(x)).(trim_ali_seqs)
    @show minimum(match_ratios), maximum(match_ratios), mean(match_ratios)
    keeps=match_ratios.>0.8
    bad_match_count=sum((!).(keeps))
    reject_seqs=vcat(reject_seqs,trim_ali_seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*"_badMatch").(nams[(!).(keeps)]))
    nams=nams[keeps]
    trim_ali_seqs=trim_ali_seqs[keeps]
    end_count = length(trim_ali_seqs)
    seq_loss = (start_count - end_count + 1) / start_count
    @show in_file, seq_loss * 100
    if length(trim_ali_seqs) > 1
        trim_ali_seqs = msa_codon_align(ref_seq, degap.(trim_ali_seqs[2:end]), scoring=score_params)
    end
    write_fasta(out_dir*in_file[1:end-6]*"_f-hmff.fasta",trim_ali_seqs,seq_names=nams)
    if length(reject_seqs) > 0
        write_fasta(out_reject_dir*in_file[1:end-6]*"_f-hmff-rejects.fasta",degap.(reject_seqs),seq_names=reject_nams)
    end
    hk_rec=[donor,visit,pool,facility,start_count,end_count-1,start_count-end_count+1,
                orf_reject_count,no_start_codon_count,no_stop_codon_count,bad_match_count]
    push!(hk,hk_rec)
end

# codon_alignment = seed_chain_align(ref=ref_seq, query=query_seq)


########################## end of function definitions ###################################


println("using Julia version: $(VERSION)")
t1 = time()

in_dir="ellpaca_collection/CAPRISA_002/env_nt_porpid/"
out_dir="ellpaca_collection/CAPRISA_002/env_nt_functionals/"
out_reject_dir="ellpaca_collection/CAPRISA_002/env_nt_functional_rejects/"
mkpath(out_dir)
mkpath(out_reject_dir)
file_names=readdir(in_dir)
file_names=file_names[(x->occursin(".fasta",x)).(file_names)]

ref_nams, ref_seqs = read_fasta("panels/consensusC_env_nt.fasta")
ref_nam=ref_nams[1]
ref_seq=ref_seqs[1]
@show ref_nam
@show ref_seq

hk=DataFrame(donor=[], visit=[], pool=[], facility=[], sequences=[], functional=[], non_functional=[],
                    orf_reject=[], no_start_codon=[], no_stop_codon=[], bad_match=[])

for file in file_names[1:end]
    pairwise_filter(in_dir,file,ref_nam,ref_seq,out_dir,out_reject_dir,hk=hk)
end

println("Totaling and writing housekeeping reports ...")
mkpath("reports/")

rec = vcat(["totals","","",""], sum.(eachcol(hk)[5:end]))
push!(hk, rec)

CSV.write(out_dir * "housekeeping_functional_filter.csv",hk)

t2 = time()
println("Functional sequences filtered in $((t2 - t1) / 60) minutes.")
