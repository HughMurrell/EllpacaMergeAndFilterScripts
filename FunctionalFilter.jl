using Pkg
Pkg.activate(".")
ENV["MPLBACKEND"] = "Agg"
using FASTX, BioSequences, StatsBase, DataFrames, CSV

########### first define some functions for later use

function degap(s)
    return replace(s,"-"=>"")
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
        

function my_write_fasta(filename, seqs;
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

#
# this next function is responsible for extracting the region from the consensus
# query that matches the reference. There are two modes for determing the reading frame:
#
# 1) normal mode: use the frame contining the longest open reading frame. return that frame
#
# 2) fallback mode: use the frame with the fewest stop codons. In this mode the matching
#    region is extracted by rearching for start and stop matches whose length and tolerance
#    is controlled by the passed parameters and these must be set depending on the molecule
#    being matched.
#
function get_matching_region_for_consensus(ref, cons; use_stop_counts=false,
                size_start=3, tol_start=0, size_stop=12, tol_stop=6)
    # translate reference
    ref_aa = String(BioSequences.translate(LongDNA{4}(degap(ref))))
    # first generate three reading frames
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
    
    @show best_reading_frame
    cons_aa = cons_aa_frames[best_reading_frame]
    cons_aa_trim = cons_aa[start:stop]
    
    if use_stop_counts
        # first choose the frame with the least stop codons
        stop_counts = (x->count('*',x)).(cons_aa_frames)
        @show stop_counts
        best_reading_frame = findmin(stop_counts)[2]
        @show best_reading_frame
        best_reading_frame = 1
        @show best_reading_frame
        cons_aa = cons_aa_frames[best_reading_frame]
    
        # now look for the start match
        start_query = ApproximateSearchQuery(LongSequence{AminoAcidAlphabet}(ref_aa[1:size_start]))
        start = findfirst(start_query, tol_start,LongSequence{AminoAcidAlphabet}(cons_aa))
        start=collect(start)[1]

        # find matching end of ref in translated consensus
        rev_ref_aa = reverse(ref_aa)
        stop_query = ApproximateSearchQuery(LongSequence{AminoAcidAlphabet}(rev_ref_aa[1:size_stop]))
        rev_cons_aa = reverse(cons_aa)
        stop = findfirst(stop_query, tol_stop, LongSequence{AminoAcidAlphabet}(rev_cons_aa))
        stop=1+length(cons_aa)-collect(stop)[1]
        cons_aa = cons_aa_frames[best_reading_frame]
        cons_aa_trim = cons_aa[start:stop]
    end

    cons_trim = cons[((start-1)*3+best_reading_frame):((stop)*3+best_reading_frame-1)]
    return cons_trim, cons_aa_trim, start, stop
end


function get_matching_region(ref, query; size=10, tol=1)
    start_query = ApproximateSearchQuery(LongDNA{4}(ref[1:size]))
    start = findfirst(start_query, tol, LongDNA{4}(query))
    # println("start = ",start)
    # stop_query = ApproximateSearchQuery(LongDNA{4}(ref[end-(size*3):end]))
    # stop = findlast(stop_query, tol*4, LongDNA{4}(query))
                                
    stop_query = ApproximateSearchQuery(LongDNA{4}(reverse(ref[end-(size):end])))
    stop = findfirst(stop_query, tol, LongDNA{4}(reverse(query)))
                                
    # println("stop = ",stop)
    if ( isnothing(start) | isnothing(stop) )
        # println("warning --- no matching regions, discarding ... ")
        return(query,translate_to_aa(query),false,"no_matching_region")
    end
    # cr = query[collect(start)[1]:collect(stop)[end]]
    cr = query[collect(start)[1]:length(query)-collect(stop)[1]+1]
    cr_aa = translate_to_aa(cr)
    if ( (length(cr) % 3) != 0 )
        return(cr,cr_aa,false,"frame_error")
    end
    # now do triplet align of trimmed query to ref and resolve
    # ali = triplet_nw_align(ref, cr, edge_reduction = 0.99, boundary_mult = 2)
    # ali = triplet_kmer_seeded_align(ref, cr,
    #    wordlength = 30, skip = 9, boundary_mult = 2,
    #    alignedcodons = true, debug=true)
    # println("ali_align = ", ali[2])
    # ali = resolve_alignments(ali[1], ali[2], mode = 1)
    # cr = degap.(ali[2])
    # println("alignment length = ",length(cr))
    # println("cr = ",cr)
    # cr = cr[1:(3*div(length(cr),3))]
    # println(cr)
    # cr_aa = String(BioSequences.translate(LongDNA{4}(cr)))
    # println("cr_aa = ",cr_aa)
    stops = findall((x->x=='*').(collect(cr_aa)))
    if length(stops) > 0 && stops[1] != length(cr_aa)
        cr = cr[1:stops[1]*3]
        cr_aa = String(BioSequences.translate(LongDNA{4}(cr)))
        # println("warning --- internal stop codon at $(stops[1]), discarding ... ")
        return(cr,cr_aa,false,"internal_stop")
    end
    frame_errs =  findall((x->x=='X').(collect(cr_aa)))
    if length(frame_errs) > 0
        # println("warning --- frame error, discarding ... ")
        return(cr,cr_aa,false,"frame_error")
    end
    if (cr_aa[1] != 'M') | (cr_aa[end] != '*')
        # println("warning --- seq has no START-and-STOP codon, discarding ... ")
        return(cr,cr_aa,false,"no_start_stop")
    end
    return(cr, cr_aa, true, "")
end



function extract_ref_match(ref_file, in_dir, file, out_dir ; house_keeping = nothing )
    
    donor=file[1:6]
    # set the output path names
    faa=mkpath(out_dir * "functional_aa/")
    out_seq_aa = faa*"/"*donor*"_functional_aa.fasta"
    faa_ali=mkpath(out_dir * "functional_aa_ali/")
    out_ali_aa = faa_ali*"/"*donor*"_functional_aa_ali.fasta"
    fnu=mkpath(out_dir * "functional_nu/")
    out_seq_nu = fnu*"/"*donor*"_functional_nu.fasta"
    fnu_ali=mkpath(out_dir * "functional_nu_ali/")
    out_ali_nu = fnu_ali*"/"*donor*"_functional_nu_ali.fasta"
    nfaa=mkpath(out_dir * "non_functional_aa/")
    nf_out_seq_aa = nfaa*"/"*donor*"_non_functional_aa.fasta"
    nfnu=mkpath(out_dir * "non_functional_nu/")
    nf_out_seq_nu = nfnu*"/"*donor*"_non_functional_nu.fasta"
    
    # first read and translate reference file
    println("processing ref and sequence ...")
    println("reading reference...")
    nr, sr = read_fasta(ref_file)
    sr_aa = String(BioSequences.translate(LongDNA{4}(degap(sr[1]))))
    if (sr_aa[1] != 'M') | (sr_aa[end] != '*')
        error("bad reference")
        exit()
    end

    my_write_fasta(out_seq_aa,[sr_aa],names=[nr[1]],aa=true,append=false)
    my_write_fasta(nf_out_seq_aa,[sr_aa],names=[nr[1]],aa=true,append=false)
    my_write_fasta(out_seq_nu,[degap(sr[1])],names=[nr[1]],aa=false,append=false)
    my_write_fasta(nf_out_seq_nu,[degap(sr[1])],names=[nr[1]],aa=false,append=false)
    println("reference read, translated and written....")
    
    # align inputfile to out_ali_nu
    println("first align the input file")
    my_mafft(in_dir * file, out_ali_nu)
    
    # get consensus of first sequence file and prepare a pannel sequence
    println("processing first visit with name $(basename(file)) ...")
    nams, seqs = read_fasta_with_names_and_descriptions(out_ali_nu)
    visits = sort(union((x->x[8:11]).(nams)))
    @show visits[1]
    cons_name="consensus of $(visits[1])"
    inds = (x->x[8:11]==visits[1]).(nams)
    ns=nams[inds]
    ss=seqs[inds]
    cons = uppercase(degap(consensus(ss)))
    println("length of reference = $(length(sr[1]))")
    println("length of consensus = $(length(cons))")
    cons, cons_aa, start, stop = get_matching_region_for_consensus(sr[1], cons)
    println("length of first visit consensus = $(length(cons))")
    println("matching region for first visit starts at $(start) and ends at $(stop)")
    
    # if start is less than 100 then try again with later visits
    for v in 2:length(visits)
        if start >= 100
            println("************************************************************************")
            @show visits[v]
            println("previous attempt no good, trying again with next visit...")
            inds = (x->x[8:11]==visits[v]).(nams)
            ns=nams[inds]
            ss=seqs[inds]
            cons = uppercase(degap(consensus(ss)))
            cons_name="consensus of $(visits[v])"
            println("length of reference = $(length(sr[1]))")
            println("length of visit $(visits[v]) consensus = $(length(cons))")
            cons, cons_aa, start, stop = get_matching_region_for_consensus(sr[1], cons)
            println("matching region for second attempt starts at $(start) and ends at $(stop)")
            println("replacing internal stops with unknowns")
            stops = findall((x->x=='*').(collect(cons_aa)))
            @show stops
            for i in 1:length(stops)-1
                cons_aa = cons_aa[1:stops[i]-1]*"-"*cons_aa[stops[i]+1:end]
                cons = cons[1:(stops[i]-1)*3]*"---"*cons[(stops[i]+1)*3-2:end]
            end
            println("************************************************************************")
        end
    end
    
    # first and second visit no good, try with first visit but use end matching
    if start >= 100
        println("************************************************************************")
        println("all visit attempts no good, trying again looking for begining and end bits...")
        inds = (x->x[8:11]==visits[1]).(nams)
        @show visits[1]
        cons_name="consensus of $(visits[1])"
        ns=nams[inds]
        ss=seqs[inds]
        cons = uppercase(degap(consensus(ss)))
        println("length of reference = $(length(sr[1]))")
        println("length of second visit consensus = $(length(cons))")
        cons, cons_aa, start, stop = get_matching_region_for_consensus(sr[1], cons, use_stop_counts=true)
        println("matching region for second attempt starts at $(start) and ends at $(stop)")
        println("replacing internal stops with unknowns")
        stops = findall((x->x=='*').(collect(cons_aa)))
        @show stops
        # cons_aa = replace(cons_aa, '\*'=>'X')
        # cons_aa = cons_aa[1:end-1]*'\*'
        for i in 1:length(stops)-1
            cons_aa = cons_aa[1:stops[i]-1]*"-"*cons_aa[stops[i]+1:end]
            cons = cons[1:(stops[i]-1)*3]*"---"*cons[(stops[i]+1)*3-2:end]
        end
        println("************************************************************************")
    end
    
    stops = findall((x->x=='*').(collect(cons_aa)))
    if length(stops) == 0 || stops[1] != length(cons_aa)
        # if its not right now, abort
        println(" ********************* still not right, aborting ***********************************")
        return()
        cons = uppercase(degap(consensus(ss)))
        cons,cons_aa=get_matching_region_for_consensus(sr[1], cons)
    end
    println("length of matching region = $(length(cons))")
    if length(cons)==0
        println("Consensus has no matching regions...")
        return()
    end
    println("consensus generated and coding region extracted...")
    my_write_fasta(out_seq_aa,[cons_aa],
        names=[cons_name], aa=true,append=true)
    # also write the nucs for later use
    my_write_fasta(out_seq_nu,[degap(cons)],
        names=[cons_name], aa=false, append=true)

    # now read all sequences and extract reading frames by aligning to cons
    # println("now extracting coding regions from seqs by comparing to consensus...")
    # println()
    for k in 1:length(visits)
        inds = (x->x[8:11]==visits[k]).(nams)
        ns=nams[inds]
        ss=seqs[inds]
        # to get rid of two 5 letter visit codes
        ns=(x->replace(x,"0a"=>"0")).(ns)
        ns=(x->replace(x,"0b"=>"1")).(ns)
        println("processing visit $(visits[k]) of length $(length(ss)) from file $(basename(file))  ")
        count_good = 0
        count_bad = 0
        count_nmr = 0
        count_isc = 0
        count_rfe = 0
        count_nss = 0
        count_length = 0
        for i in 1:length(ss)
            ts = uppercase(degap(ss[i]))
            ts, ts_aa, ok, msg = get_matching_region(cons, ts, size=15, tol=2)
            if ok
                count_length += length(ts_aa)
                my_write_fasta(out_seq_aa,[ts_aa],
                    names=[ns[i]], aa=true, append=true)
                my_write_fasta(out_seq_nu,[degap(ts)],
                    names=[ns[i]], aa=false, append=true)
                count_good += 1
            else
                my_write_fasta(nf_out_seq_aa,[ts_aa],
                    names=[ns[i]*" "*msg], aa=true, append=true)
                my_write_fasta(nf_out_seq_nu,[degap(ts)],
                    names=[ns[i]*" "*msg], aa=false, append=true)
                count_bad += 1
                msg == "no_matching_region" ? count_nmr += 1 : nothing
                msg == "internal_stop" ? count_isc += 1 : nothing
                msg == "frame_error" ? count_rfe += 1 : nothing
                msg == "no_start_stop" ? count_nss += 1 : nothing
            end
        end
        print("processed $(count_good) good seqs and $(count_bad) bad seqs, ")
        println("success rate = $(count_good / (count_good + count_bad))")
        if ! isnothing(house_keeping)
            visit=basename(file)[8:11]
            mean_aa_length = count_good > 0 ? Int(round(count_length / count_good)) : 0
            push!(house_keeping,[donor,visits[k],
                        count_good,count_bad,
                        mean_aa_length,count_nmr,count_isc,
                        count_rfe,count_nss])
        end
    end
    println()
    println("using mafft to align aa sequences (please be patient)....")
    my_mafft(out_seq_aa, out_ali_aa)
    
    println("generating codon aware nucleotide alignment ...")
    nr, sr = read_fasta_with_names_and_descriptions(out_ali_aa)
    nu_n, nu_s = read_fasta_with_names_and_descriptions(out_seq_nu)
    
    rfe_count = 0
    rm(out_ali_nu,force=true)
    for i in 1:length(nu_s)
        re_gapped_ns = reverse_translate_with_gaps(sr[i]*"*",nu_s[i])
        if length(re_gapped_ns) > 0
            my_write_fasta(out_ali_nu,[re_gapped_ns], names=[nu_n[i]], aa=false, append=true)
        else
            println(" ****** $(nu_n[i]) has posthumus reading frame error, this should not happen ******")
            # println("appending $(nu_n[i]) to $(nf_out_seq_nu)")
            rfe_count += 1
            my_write_fasta(nf_out_seq_nu,[nu_s[i]], names=[nu_n[i] * " reading_frame_error"], append=true)
        end
    end
    house_keeping[end,:reading_frame_error] += rfe_count
    house_keeping[end,:non_functional] += rfe_count
    house_keeping[end,:functional] -= rfe_count
    
    println("sorting rejects on reason .....")
    nu_n, nu_s = read_fasta_with_names_and_descriptions(nf_out_seq_nu)
    my_write_fasta(nf_out_seq_nu,[nu_s[1]], names=[nu_n[1]], aa=false, append=false)
    if length(nu_n) > 1
        nu_n=nu_n[2:end]
        nu_s=nu_s[2:end]
        reasons=(s->split(s," ")[4]).(nu_n)
        p=sortperm(reasons)
        nu_n=nu_n[p]
        nu_s=nu_s[p]
        my_write_fasta(nf_out_seq_nu,nu_s, names=nu_n, aa=false, append=true)
    end
    aa_n, aa_s = read_fasta_with_names_and_descriptions(nf_out_seq_aa)
    my_write_fasta(nf_out_seq_aa,[aa_s[1]], names=[aa_n[1]], aa=true, append=false)
    if length(aa_n) > 1
        aa_n=aa_n[2:end]
        aa_s=aa_s[2:end]
        reasons=(s->split(s," ")[4]).(aa_n)
        p=sortperm(reasons)
        aa_n=aa_n[p]
        aa_s=aa_s[p]
        my_write_fasta(nf_out_seq_aa,aa_s, names=aa_n, aa=true, append=true)
    end
    println("$(donor) completed ...")
end

########################## end of function definitions ###################################


println("using Julia version: $(VERSION)")
t1 = time()

in_dir="ellpaca_nu_merged/"
out_dir="ellpaca_functionals/"
mkpath(out_dir)
file_names=readdir(in_dir)
file_names=file_names[(x->occursin("nu_merged.fasta",x)).(file_names)]
donors=sort((x->x[1:6]).(file_names))
@show donors
println("Number of donors = $(length(donors))")

fr="data_in/panels/hxb2-env.fasta"
println(" doing codon aware alignment against reference in $(fr) ...")

hk=DataFrame(donor=[], visit=[], functional=[],non_functional=[],
                        mean_functional_aa_length=[],
                        no_matching_region=[],internal_stop_codon=[],
                        reading_frame_error=[],no_stop_and_start=[])
                        
for file in file_names[1:end]
    @show(basename(file))
    extract_ref_match(fr, in_dir, file, out_dir, house_keeping=hk)
    println()
end

println("Totaling and writing housekeeping reports ...")
mkpath("reports/")

ghk = groupby(hk, [:donor])
sghk = combine(ghk, [:functional => sum, :non_functional => sum,
                     :no_matching_region => sum, :internal_stop_codon => sum,
                     :reading_frame_error => sum, :no_stop_and_start => sum ], renamecols=false)
                     
rec = vcat(["totals",""], sum.(eachcol(hk)[3:end]))
rec[5]=mode(hk[!,:mean_functional_aa_length])
push!(hk, rec)

CSV.write("reports/functional_filter_visit_housekeeping.csv",hk)
                     
rec = vcat(["totals"], sum.(eachcol(sghk)[2:end]))
push!(sghk, rec)
                     
CSV.write("reports/functional_filter_donor_housekeeping.csv",sghk)

t2 = time()
println("Functional sequences filtered in $((t2 - t1) / 3600) hours.")
