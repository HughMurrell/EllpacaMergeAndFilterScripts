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

function rename_file_and_seqs(in_dir, out_dir, in_file_name, chk)
    nams, seqs = read_fasta_with_names_and_descriptions(in_dir*in_file_name)
    donor=basename(in_file_name)[1:6]
    visit=basename(in_file_name)[8:11]
    # do some housekeeping
    colnames = names(chk)
    visits=union((x->split(x,"_")[2][1:4]).(nams))
    for visit in visits
        visit_nams = (x->(split(x,"_")[2][1:4]==visit)).(nams)
        pool_facility="p"*split(in_file_name,"-")[1][5:end]
        rec=[donor,visit,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        for i in 1:length(colnames)
            if colnames[i] == pool_facility
                rec[i]=sum(visit_nams)
            end
        end
        push!(chk,rec)
    end
    pool=split(in_file_name,"_")[1][5:end]
    while length(pool)<2
        pool="0"*pool
    end
    facility=split(in_file_name,"_")[2][1:4]
    facility=split(facility,"-")[1]
    out_file_name = "C002_"*donor*"_"*visit*"_env_nt_s-0000-p"*pool*"-"*facility*"_p-fs5-ma0p7"
    if donor*"_"*visit in SUSPECTS
        out_file_name=out_file_name * "_e-suspect"
    end
    out_file_name = out_file_name*".fasta"
    new_nams = (x -> out_file_name[1:20] * "_pb-" * x[12:19] * "_p-fs" *
            split(split(x," ")[2],"=")[2] * "-ma" * split(split(x," ")[3],"=")[2] ).(nams)
    write_fasta(out_dir*out_file_name,degap.(seqs),names=new_nams)
    @show out_file_name, length(seqs)
end


############## end of function definitions

println("Rename script")
println("using Julia version: $(VERSION)")
t1 = time()
pools=[1,2,3,4,5,6,7,8,8,9,9,10,10,11,12]
# pools=[8,8]
in_dir="data_in/"
out_dir="ellpaca_collection/CAPRISA_002/env_nt_porpid/"
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
        
SUSPECTS = [ "CAP045_2060",
             "CAP210_2050",
             "CAP265_4170",
             "CAP283_4220",
             "CAP296_4240",
             "CAP314_4190",
             "CAP322_6110",
             "CAP327_4180",
             "CAP382_3090"]

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

for pool in pool_dirs
    for sample in pool_sample_dict[pool]
        fs = pool*"/"*sample*"/"*sample*".fasta"
        rename_file_and_seqs(in_dir,out_dir,fs,c_hk)
    end
end

# sort the housekeeping table before totalling
sort!(c_hk, [:donor, :visit])

# tally up the file counts:
out_files=readdir(out_dir)
out_files=out_files[(s -> endswith(s,".fasta")).(out_files)]
out_files=(x->join(split(x,"_")[2:3],"_")).(out_files)
hk_files=c_hk[:,:donor].*"_".*c_hk[:,:visit]
@show length(out_files), length(hk_files)
@show setdiff(out_files,hk_files)
@show setdiff(hk_files,out_files)

cm = countmap(hk_files)
for k in keys(cm)
    if cm[k] > 2
        @show k, cm[k]
    end
end

# now sum up and save the housekeeping tables

rec = vcat(["totals","all"], sum.(eachcol(c_hk)[3:end]))
push!(c_hk, rec)
tot_col = sum.(eachrow(c_hk[:, names(c_hk, Real)]))
c_hk[!,"totals"]=tot_col
mkpath("reports/")
CSV.write(out_dir * "porpid_housekeeping.csv",c_hk)

t2 = time()
println("All Renaming Done in $((t2 - t1)/60) minutes, see porpid_housekeeping.csv")
exit()

