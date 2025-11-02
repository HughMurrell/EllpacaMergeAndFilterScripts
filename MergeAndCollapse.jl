using Pkg
Pkg.activate(".")
ENV["MPLBACKEND"] = "Agg"
using FASTX, BioSequences, StatsBase, DataFrames, CSV
using MAFFT_jll

########### first define some functions for later use

function degap(s::String)
    return replace(s,"-"=>"")
end

function mafft(inpath1, inpath2, outpath)
    print("profile aligning ....")
    cmd = `$(mafft_profile()) --quiet --thread 2 --ep 2 --op 3 --out $(outpath) $(inpath1) $(inpat2)`
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

function merge_seqs_by_donor(in_dir, out_dir, donor, files)

    # first we collect
    write_fasta(out_dir * donor * "_nt_collected.fasta",[], names=[],append=false)

    for file in files
        nams, seqs = read_fasta_with_names_and_descriptions(in_dir * file)
        # remove ref and degap
        keeps=(x->startswith(x,"C002")).(nams)
        nams=nams[keeps]
        seqs=degap.(seqs[keeps])
        println("$(file) -> $(length(seqs)) sequences")
        write_fasta(out_dir * donor * "_nt_collected.fasta", seqs, names=nams, append=true)
    end
    
    # now we merge
    nams, seqs = read_fasta_with_names_and_descriptions(out_dir * donor * "_nt_collected.fasta")
    seqs=(degap).(seqs)
    visits=sort(union((x->split(x,"_")[3]).(nams)))
    
    # sort nams and seqs prior to merge
    sp=sortperm(nams)
    nams=nams[sp]
    seqs=seqs[sp]
    
    prim_nams=(x->join(split(x,"_")[1:5],"_")).(nams)
    # extra_nams=(x->join(split(x,"_")[7:end],"_")).(nams)
    
    p_nams=(x->split(x,"_")[6]).(nams)
    fs_ints=(x->parse(Int, String(split(x,"-")[2])[3:end])).(p_nams)
    ma_floats=(x->parse(Float64, String(split(x,"-")[3])[3:end])).(p_nams)
    
    # dropped the -merg from the f flag on Chris request
    out_file=out_dir * "C002_" * donor * "_" * join(visits,"-") * "_env_nt_p-fs5-ma0p7_f-hmff.fasta"
    write_fasta(out_file,[], names=[],append=false)
    
    r=1
    n=2
    while n <= length(nams)
        if ( prim_nams[n] == prim_nams[r] ) && ( seqs[n] != seqs[r] )
            println(prim_nams[n]," has dif seq to ",prim_nams[r])
            if fs_ints[n] > fs_ints[r]
                println("************ using better sequence ******************************")
                seqs[r] = seqs[n]
                fs_ints[r] = fs_ints[n]
                println("*****************************************************************")
            end
            n=n+1
        end
        if ( prim_nams[n] == prim_nams[r] ) && ( seqs[n] == seqs[r] )
            fs_ints[r] = fs_ints[r]+fs_ints[n]
            ma_floats[r] = min(ma_floats[r],ma_floats[n])
            n=n+1
        else
            write_fasta(out_file,[seqs[r]],
                names=[prim_nams[r]*"_p-fs$(fs_ints[r])-ma$(ma_floats[r])"],
                append=true)
            # global seq_out_count += 1
            n=n+1
            r=n-1
        end
    end
    while r<=length(nams)
        write_fasta(out_file,[seqs[r]],
                names=[prim_nams[r]*"_p-fs$(fs_ints[r])-ma$(ma_floats[r])"],
                append=true)
        # global seq_out_count += 1
        r=r+1
    end
    
    # now remove the collected file
    rm(out_dir * donor * "_nt_collected.fasta")
    
end

function variant_collapse(seqs; prefix = "seq_")
    dic = countmap(uppercase.(seqs))
    OC = sort(dic, by = x -> dic[x], rev = true)
    seqs,sizes = collect(keys(OC)),collect(values(OC))
    names = [prefix*string(i)*"-"*string(sizes[i]) for i in 1:length(seqs)]
    return seqs,sizes,names
end


############## end of function definitions ##############################

println("MERGE script")
println("using Julia version: $(VERSION)")
t1 = time()
in_dir="ellpaca_collection/CAPRISA_002/env_nt_functionals/"
out_dir="ellpaca_collection/CAPRISA_002/env_nt_merged/"
mkpath(out_dir)

#now lets see what we have....

donor_files=readdir(in_dir)
donor_files=donor_files[(s -> endswith(s,".fasta")).(donor_files)]
donor_files=donor_files[(s -> ! occursin("suspect",s)).(donor_files)]
donors=sort(union((x->split(x,"_")[2]).(donor_files)))
println("Number of donors = $(length(donors))")

for donor in donors
    fs = []
    for file in donor_files
        if split(file,"_")[2] == donor
            push!(fs,file)
        end
    end
    println()
    si=(x->split(x,"_")[3]).(fs)
    fs=fs[sortperm(si)]
    println("------- $(donor) -------")
    merge_seqs_by_donor(in_dir, out_dir, donor, fs)
end

# now some housekeeping for the merged sequences
file_names=readdir(out_dir)
file_names=file_names[(x->occursin("merg",x)).(file_names)]
all_visits=[]

for file in file_names
    nams, seqs = read_fasta_with_names_and_descriptions(out_dir * file)
    visits=union((x->split(x,"_")[3][1:4]).(nams))
    global all_visits=union(vcat(all_visits,visits))
end

sort!(all_visits)

println("doing merge housekeeping ....")

m_hk = DataFrame( [ [] for i in 1:length(all_visits)+1 ],
                vcat( ["donor"], all_visits ) )

for file in file_names
    nams, seqs = read_fasta_with_names_and_descriptions(out_dir * file)
    visits=union((x->split(x,"_")[3][1:4]).(nams))
    donor=split(file,"_")[2]
    next_rec=vcat([donor],zeros(Int,length(all_visits)))
    push!(m_hk,next_rec)
    for visit in visits
        hits=sum( (x->split(x,"_")[3][1:4]).(nams) .== visit )
        m_hk[end,visit]=hits
    end
end
rec = vcat(["totals"], sum.(eachcol(m_hk)[2:end]))
push!(m_hk, rec)

# final housekeepig
cp_hk = DataFrame( donor = String[], merged = Int[], collapsed = Int[] )

println("generating collapsed merged sequences...")
in_dir = out_dir
out_dir = "ellpaca_collection/CAPRISA_002/env_nt_collapsed/"
mkpath(out_dir)
file_names=readdir(in_dir)
file_names=file_names[(x->occursin(".fasta",x)).(file_names)]
merged_tot=0
collapsed_tot=0
for file in file_names
    in_file = in_dir * file
    out_file = out_dir * file[1:end-6] * "-coll.fasta"
    println(" collapsing $(in_file)" )
    nams, seqs = read_fasta_with_names_and_descriptions(in_file)
    merged_count = length(seqs)
    collapsed_count = 0
    seqs=(degap).(seqs)
    visits=sort(union((x->split(x,"_")[3][1:4]).(nams)))
    write_fasta(out_file,[], names=[],append=false)
    for visit in visits
        visit_inds = (x->(split(x,"_")[3][1:4]==visit)).(nams)
        visit_nams = nams[visit_inds]
        visit_seqs = seqs[visit_inds]
        # collapse each visit segment
        pref = "C002_" * split(file,"_")[2] * "_" * visit * "_env_pb-coll-"
        col_seqs, col_sizes, col_names = variant_collapse(visit_seqs, prefix = pref)
        write_fasta(out_file, col_seqs, names = col_names, append=true)
        collapsed_count += length(col_seqs)
    end
    next_rec = ["$(split(file,"_")[2])", merged_count, collapsed_count]
    push!(cp_hk, next_rec)
    global merged_tot += merged_count
    global collapsed_tot += collapsed_count
end
rec = ["totals", merged_tot, collapsed_tot]
push!(cp_hk, rec)

mcp_hk = innerjoin(m_hk, cp_hk, on = :donor)
CSV.write(out_dir * "housekeeping_merged_collapsed.csv",mcp_hk)

t2 = time()
println()
println("All functional files merged and collapsed in $((t2 - t1)/60) minutes.")

