using SeisIO, SeisNoise, JLD2, Distributed

@everywhere include("../../src/reference.jl")
include("../../src/convergence.jl")
include("../../src/utils.jl")

# project directory
base_dir = "/Users/jared/SCECintern2019/data"
dir = "BPnetwork_2004"

# input cross-correlations
corrname  = "$base_dir/xcorrs/$dir/BPnetwork_2004_xcorrs"
# reference output file
base_dir = "/Users/jared/SCECintern2019/data"
# input references
refname = "$base_dir/reference/$dir/BPnetwork_2004_ref.jld2"

# input file holds metadata (stationlist, timestamplist, etc.)
f = jldopen(corrname*".jld2")
tslist = f["info/timestamplist"]
close(f) # base xcorr file

# generate reference by stacking all cross-correlations
ref_dicts = pmap(t->compute_reference_xcorr(t, corrname, phase_smoothing=0., stack="linear"), tslist)

# collect all references into one dictionary
ref_dict_out = Dict()
for i=1:length(ref_dicts)
    for stnpair in keys(ref_dicts[i])
        if haskey(ref_dict_out, stnpair)
            ref_dict_out[stnpair].corr .+= ref_dicts[i][stnpair].corr
        else
            ref_dict_out[stnpair] = copy(ref_dicts[i][stnpair])
        end
    end
end

f_out = jldopen(refname, "w")
for stnpair in keys(ref_dict_out)
    f_out[stnpair] = ref_dict_out[stnpair]
end
close(f_out)
rmprocs(workers())

# compute convergence of successively long stacks to the reference xcorr
# coda wave
convname_cc = "$base_dir/convergence/$dir/BPnet_2004_coda_sel.jld2"
convergence(corrname, refname2, convname_cc, ntimes=90, metric="snr", slice=[10.0, 50.0], thresh=0.0)
# ballistic wave
convname_cc = "$base_dir/convergence/$dir/BPnet_2004_bal_sel.jld2"
convergence(corrname, refname2, convname_cc, ntimes=90, metric="snr", slice=10.0, thresh=0.0)
