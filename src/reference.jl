include("utils.jl")
include("io.jl")
include("stacking.jl")

export compute_reference_xcorr

function compute_reference_xcorr()
    # computing reference xcorr
    # input cross-correlations
    corrname  = Outputdir*"/cc/$(Year)_xcorrs"
    # output references
    refname = corrname*"_ref.jld2" # this is fixed in the SeisXcorrelation/pstack.

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

end

"""
    map_reference(tstamp::String, basefiname::String; phase_smoothing::Float64=0., stack::String="linear", reference::String="", thresh::Float64=-1.0)

Stack all cross-correlation functions for all given station pairs to generate a reference cross-correlation.

# Arguments
- `tstamp::String`    :
- `basefiname::String,`    : Input base file name e.g. "./inputData/BPnetwork"
- `phase_smoothing::Float64`    : Level of phase_smoothing (0 for linear stacking)
- `stack::String`     : "selective" for selective stacking and "linear" for linear stacking
- `reference::String`    : Path to the reference used in selective stacking
- `thresh::Float64`     : Threshold used for selective stacking

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair
"""
function map_reference(tstamp::String, basefiname::String; phase_smoothing::Float64=0., stack::String="linear", reference::String="", thresh::Float64=-1.0)
    # hold reference xcorrs in memory and write all at once
    ref_dict = Dict()

    # read unstacked xcorrs for each time stamp
    f_cur = jldopen(basefiname*".$tstamp.jld2")
    grp = f_cur[tstamp] # xcorrs
    println("$tstamp")

    # iterate over station pairs
    for pair in sort(keys(grp))
        # TODO: use only unique station pairings when creating references. Currently no guarantee of uniqueness (reverse can exist)
        # load xcorr
        xcorr = try grp[pair] catch; continue end
        #remove_nan!(xcorr)

        # stack xcorrs over length of CorrData object using either "selective" stacking or "linear" stacking
        if stack=="selective"
            f_exref = jldopen(reference)
            ref = f_exref[pair]
            close(f_exref)
            xcorr, rmList = selective_stacking(xcorr, ref, threshold=thresh)
        else
            stack!(xcorr, allstack=true, phase_smoothing=phase_smoothing)
        end

        # stack xcorrs if they have a key, assign key if not
        if haskey(ref_dict, pair)
            ref_dict[pair].corr .+= xcorr.corr
        else
            ref_dict[pair] = xcorr
        end
    end
    close(f_cur) # current xcorr file

    return ref_dict
end
