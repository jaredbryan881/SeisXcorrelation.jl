using SeisIO, JLD2, SeisNoise, Statistics, PlotlyJS
include("utils.jl")
include("io.jl")
include("stacking.jl")

"""
    compute_reference_xcorr(tstamp::String, basefiname::String; phase_smoothing::Float64=0., stack::String="linear", reference::String="", thresh::Float64=-1.0)

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
function compute_reference_xcorr(tstamp::String, basefiname::String; phase_smoothing::Float64=0., stack::String="linear", reference::String="", thresh::Float64=-1.0)
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
