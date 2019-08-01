using SeisIO, JLD2, SeisNoise, Statistics, PlotlyJS
include("utils.jl")
include("io.jl")
include("stacking.jl")

"""
    compute_reference_xcorr(finame::String, foname::String, phase_smoothing::Float64=0., stack::String="selective")

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `basefiname::String,`    : Input base file name e.g. "./inputData/BPnetwork"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"
- `phase_smoothing::Float64`    : Level of phase_smoothing (0 for linear stacking)
- `stack::String`     : "selective" for selective stacking and "linear" for linear stacking

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair
"""
function compute_reference_xcorr(basefiname::String, foname::String; phase_smoothing::Float64=0., stack::String="selective")
    # input file holds metadata (stationlist, timestamplist, etc.)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"]
    close(f) # base xcorr file

    # hold reference xcorrs in memory and write all at once
    ref_dict = Dict()

    # iterate over timestamps
    for (t, tstamp) in enumerate(tslist)
        # read unstacked xcorrs for each time stamp
        f_cur = jldopen(basefiname*".$tstamp.jld2")
        grp = f_cur[tstamp] # xcorrs
        println("stacking xcorrs at $tstamp")

        # iterate over station pairs
        for pair in keys(grp)
            # load xcorr and reverse if necessary
            xcorr = try grp[pair] catch; continue end

            # stack xcorrs over length of CorrData object using either "selective" stacking or "linear" stacking
            if stack=="selective"
                xcorr, rmList = selective_stacking(xcorr)
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
    end

    # write output file of reference xcorrs
    f_ref = jldopen(foname, "a+")
    for pair in keys(ref_dict)
        f_ref[pair] = ref_dict[pair]
    end
    close(f_ref) # reference xcorr file
end

"""
    compute_reference_xcorr(finame::String, existing_reference::String, foname::String; phase_smoothing::Float64=0.)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `basefiname::String,`    : Input base file name e.g. "./inputData/BPnetwork"
- `existing_reference::String`    : Input file name for precomputed reference to improve e.g. "./refData/referenceXcorr.jld2"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"
- `phase_smoothing::Float64=0.`    : Level of phase smoothing for phase weighted stacking. 0. for linear stacking.

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair
"""
function compute_reference_xcorr(basefiname::String, existing_reference::String, foname::String; phase_smoothing::Float64=0., threshold::Float64=0.0)
    # input file holds metadata (stationlist, timestamplist, etc.)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"]
    close(f) # base xcorr file

    # hold reference xcorrs in memory and write all at once
    ref_dict = Dict()

    # existing reference to compare windows to
    f_ref = jldopen(existing_reference)

    # iterate over timestamps
    for (t, tstamp) in enumerate(tslist)
        # read unstacked xcorrs for each time stamp
        f_cur = jldopen(basefiname*".$tstamp.jld2")
        grp = f_cur[tstamp] # xcorrs
        println("stacking xcorrs at $tstamp")

        # iterate over station pairs
        for pair in keys(grp)
            # load xcorr
            xcorr = try grp[pair] catch; continue end
            # load reference xcorr
            ref = try f_ref[pair] catch; continue end

            # stack xcorrs over length of CorrData object using selective stacking
            xcorr, nRemoved = selective_stacking(xcorr, ref, threshold=threshold)

            # stack xcorrs if they have a key, assign key if not
            if haskey(ref_dict, pair)
                ref_dict[pair].corr[:,1] .+= xcorr.corr[:,1]
            else
                ref_dict[pair] = xcorr
            end
        end
        close(f_cur) # current xcorr file
    end

    # write output file of reference xcorrs
    f_ref = jldopen(foname, "a+")
    for pair in keys(ref_dict)
        f_ref[pair] = ref_dict[pair]
    end
    close(f_ref) # reference xcorr file
end
