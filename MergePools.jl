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

function collect_and_merge_seqs_by_donor(in_dir, out_dir, donor, files, chk)
    
    colnames = names(chk)
    
    # first we collect
    write_fasta(out_dir * donor * "_nu_collected.fasta",[], names=[],append=false)

    for file in files
        nams, seqs = read_fasta_with_names_and_descriptions(in_dir * file)
        println("$(file) -> $(length(seqs)) sequences")
        write_fasta(out_dir * donor * "_nu_collected.fasta", seqs, names=nams, append=true)
        
        # do some housekeeping
        visits=union((x->split(x,"_")[2][1:4]).(nams))
        for visit in visits
            visit_nams = (x->(split(x,"_")[2][1:4]==visit)).(nams)
            pool_facility="p"*split(file,"-")[1][5:end]
            rec=[donor,visit,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            for i in 1:length(colnames)
                if colnames[i] == pool_facility
                    rec[i]=sum(visit_nams)
                end
            end
            push!(chk,rec)
        end
    end
    
    # now we merge
    nams, seqs = read_fasta_with_names_and_descriptions(out_dir * donor * "_nu_collected.fasta")
    seqs=(degap).(seqs)
    visits=union((x->split(x,"_")[2][1:4]).(nams))
    
    # sort nams and seqs prior to merge
    sp=sortperm(nams)
    nams=nams[sp]
    seqs=seqs[sp]
    
    prim_nams=(x->split(x," ")[1]).(nams)
    
    ccs_nams=(x->split(x," ")[2]).(nams)
    ccs_ints=(x->parse(Int, String(split(x,"=")[2]))).(ccs_nams)
    ma_nams=(x->split(x," ")[3]).(nams)
    ma_floats=(x->parse(Float64, String(split(x,"=")[2]))).(ma_nams)
    
    out_file=out_dir * donor * "_nu_merged.fasta"
    write_fasta(out_file,[], names=[],append=false)
    
    r=1
    n=2
    while n <= length(nams)
        if ( prim_nams[n] == prim_nams[r] ) && ( seqs[n] != seqs[r] )
            println(prim_nams[n]," has dif seq to ",prim_nams[r])
            if ccs_ints[n] > ccs_ints[r]
                println("************ using better sequence ******************************")
                seqs[r] = seqs[n]
                ccs_ints[r] = ccs_ints[n]
                println("*****************************************************************")
            end
            n=n+1
        end
        if ( prim_nams[n] == prim_nams[r] ) && ( seqs[n] == seqs[r] )
            ccs_ints[r] = ccs_ints[r]+ccs_ints[n]
            ma_floats[r] = min(ma_floats[r],ma_floats[n])
            n=n+1
        else
            write_fasta(out_file,[seqs[r]],
                names=[prim_nams[r]*" fs=$(ccs_ints[r]) ma=$(ma_floats[r])"],
                append=true)
            # global seq_out_count += 1
            n=n+1
            r=n-1
        end
    end
    while r<=length(nams)
        write_fasta(out_file,[seqs[r]],
                names=[prim_nams[r]*" fs=$(ccs_ints[r]) ma=$(ma_floats[r])"],
                append=true)
        # global seq_out_count += 1
        r=r+1
    end
    
    # now remove the collected file
    rm(out_dir * donor * "_nu_collected.fasta")
    
end

function variant_collapse(seqs; prefix = "seq_")
    dic = countmap(uppercase.(seqs))
    OC = sort(dic, by = x -> dic[x], rev = true)
    seqs,sizes = collect(keys(OC)),collect(values(OC))
    names = [prefix*string(i)*"_"*string(sizes[i]) for i in 1:length(seqs)]
    return seqs,sizes,names
end


############## end of function definitions

println("MERGE POOL script")
println("using Julia version: $(VERSION)")
t1 = time()
pools=[1,2,3,4,5,6,7,8,8,9,9,10,10,11,12]
# pools=[8,8]
in_dir="data_in/"
out_dir="ellpaca_nu_merged/"
mkpath(out_dir)
pool_dirs=["pool$(i)_" for i in pools]
pool_dirs_ext=["nicd-postproc",
           "nicd-postproc",
           "earlham-postproc",
           "earlham-postproc",
           "uw-postproc",
           "nicd-postproc",
           "nicd-postproc",
           "nicd-postproc",
           "inqaba-postproc",
           "uw-postproc",
           "inqaba-postproc",
           "uw-postproc",
           "inqaba-postproc",
           "inqaba-postproc",
           "inqaba-postproc"]
pool_dirs=pool_dirs .* pool_dirs_ext
println(pool_dirs)

# collection housekeeping
colnames=["p1_nicd","p2_nicd","p3_earlham","p4_earlham",
        "p5_uw","p6_nicd","p7_nicd","p8_nicd",
        "p8_inqaba","p9_uw","p9_inqaba","p10_uw","p10_inqaba","p11_inqaba","p12_inqaba"]
c_hk=DataFrame(donor=String[],visit=String[],p1_nicd=Int[],p2_nicd=Int[],p3_earlham=Int[],
        p4_earlham=Int[],p5_uw=Int[],p6_nicd=Int[],p7_nicd=Int[],p8_nicd=Int[],p8_inqaba=Int[],
        p9_uw=Int[],p9_inqaba=Int[],p10_uw=Int[],p10_inqaba=Int[],p11_inqaba=Int[],p12_inqaba=Int[])

#now lets see what we have....

pool_sample_dict=Dict()
for d in pool_dirs
    samples = readdir(in_dir*d)
    samples = samples[(s -> startswith(s,"CAP")).(samples)]
    println("$(d) ======> $(length(samples)) samples")
    pool_sample_dict[d]=samples
end

donors=[]
for k in keys(pool_sample_dict)
    # println("$(k) -> $(length(pool_sample_dict[k])) samples")
    pd=(s -> s[1:6]).(pool_sample_dict[k])
    global donors = vcat(donors,pd)
end


donors=sort(union(donors))
println("Number of donors = $(length(donors))")

for donor in donors
    fs = []
    for pool in pool_dirs
        for sample in pool_sample_dict[pool]
            if sample[1:6]==donor
                push!(fs,pool*"/"*sample*"/"*sample*".fasta")
            end
        end
    end
    println()
    si=(x->split(x,"/")[3]).(fs)
    fs=fs[sortperm(si)]
    println("------- $(donor) -------")
    mhk=nothing
    collect_and_merge_seqs_by_donor(in_dir, out_dir, donor, fs, c_hk)
end

# now sum up and save the housekeeping tables

rec = vcat(["totals","all"], sum.(eachcol(c_hk)[3:end]))
push!(c_hk, rec)
tot_col = sum.(eachrow(c_hk[:, names(c_hk, Real)]))
c_hk[!,"totals"]=tot_col
mkpath("reports/")
CSV.write("reports/donor_collection_housekeeping.csv",c_hk)

# now some housekeeping for the merged sequences
file_names=readdir(out_dir)
file_names=file_names[(x->occursin("merged",x)).(file_names)]
all_visits=[]

for file in file_names
    nams, seqs = read_fasta_with_names_and_descriptions(out_dir * file)
    visits=union((x->split(x,"_")[2][1:4]).(nams))
    global all_visits=union(vcat(all_visits,visits))
end

sort!(all_visits)

println("doing merge housekeeping ....")

m_hk = DataFrame( [ [] for i in 1:length(all_visits)+1 ],
                vcat( ["donor"], all_visits ) )

for file in file_names
    nams, seqs = read_fasta_with_names_and_descriptions(out_dir * file)
    visits=union((x->split(x,"_")[2][1:4]).(nams))
    donor=file[1:6]
    next_rec=vcat([donor],zeros(Int,length(all_visits)))
    push!(m_hk,next_rec)
    for visit in visits
        hits=sum( (x->split(x,"_")[2][1:4]).(nams) .== visit )
        m_hk[end,visit]=hits
    end
end
rec = vcat(["totals"], sum.(eachcol(m_hk)[2:end]))
push!(m_hk, rec)

# final housekeepig
cp_hk = DataFrame( donor = String[], merged = Int[], collapsed = Int[] )

println("generating collapsed merged sequences...")
in_dir = out_dir
out_dir = "ellpaca_nu_merged_collapsed/"
mkpath(out_dir)
file_names=readdir(in_dir)
file_names=file_names[(x->occursin("nu_merged.fasta",x)).(file_names)]
merged_tot=0
collapsed_tot=0
for file in file_names
    in_file = in_dir * file
    out_file = out_dir * file[1:end-6] * "_collapsed.fasta"
    println(" collapsing $(in_file)" )
    nams, seqs = read_fasta_with_names_and_descriptions(in_file)
    merged_count = length(seqs)
    collapsed_count = 0
    seqs=(degap).(seqs)
    visits=sort(union((x->split(x,"_")[2][1:4]).(nams)))
    write_fasta(out_file,[], names=[],append=false)
    for visit in visits
        visit_inds = (x->(split(x,"_")[2][1:4]==visit)).(nams)
        visit_nams = nams[visit_inds]
        visit_seqs = seqs[visit_inds]
        # collapse each visit segment
        pref = file[1:6] * "_" * visit * "_v"
        col_seqs, col_sizes, col_names = variant_collapse(visit_seqs, prefix = pref)
        write_fasta(out_file, col_seqs, names = col_names, append=true)
        collapsed_count += length(col_seqs)
    end
    next_rec = ["$(file[1:6])", merged_count, collapsed_count]
    push!(cp_hk, next_rec)
    global merged_tot += merged_count
    global collapsed_tot += collapsed_count
end
rec = ["totals", merged_tot, collapsed_tot]
push!(cp_hk, rec)

mcp_hk = innerjoin(m_hk, cp_hk, on = :donor)
CSV.write("reports/donor_merged_collapsed_housekeeping.csv",mcp_hk)

t2 = time()
println()
println("All pools merged and collapsed in $((t2 - t1)/60) minutes.")

