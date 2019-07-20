using SeisIO, JLD2, SeisNoise, PlotlyJS
include("utils.jl")
include("io.jl")

"""
compute_reference_xcorr(finame::String, foname::String)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `basefiname::String,`    : Input base file name e.g. "./inputData/BPnetwork"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair

"""
function compute_reference_xcorr(basefiname::String, foname::String)
    # input file holds metadata (stationlist, timestamplist, etc.)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"]
    close(f) # base xcorr file

    # hold reference xcorrs in memory and write all at once
    ref_dict = Dict()

    # iterate over timestamps
    for tstamp in tslist
        # read unstacked xcorrs for each time stamp
        f_cur = jldopen(basefiname*".$tstamp.jld2")
        grp = f_cur[tstamp] # xcorrs
        println("stacking xcorrs at $tstamp")

        # iterate over station pairs
        for pair in keys(grp)
            # load xcorr and reverse if necessary
            xcorr = grp[pair]

            # stack xcorrs over DL_time_unit (or over time length of CorrData object)
            stack!(xcorr, allstack=true)

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
convergence(finame::String, foname::String; metric::String="snr", skipoverlap::Bool=false)

Compute the RMS error or SNR for progressively longer stacks of cross-correlations and (for RMS) a reference cross-correlation

# Arguments
- `finame::String,`    : Input file name for unstacked cross-correlations e.g. "./inputData/BPnetwork.jld2"
- `refname::String,`    : Input file name for reference cross-corrleations e.g. "./inputData/BPnetwork.jld2"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"
- `metric::String`    : Method of tracking convergence: "snr" for peak signal-to-noise ratio or "rms" for root mean square misfit
- `skipoverlap::Bool`    : Whether or not to skip overlap in cross-correlation windows. If cross-correlations were computed
                           with no overlap, set this to false to ensure no time is skipped.

# Output
- `foname.jld2`    : contains arrays of RMS errors for progressively longer stacks of cross-correlations

"""
function convergence(basefiname::String, refname::String, foname::String; metric::String="snr", skipoverlap::Bool=false)
    # input file contains metadata (stationlist, timestamplist)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"] # time stamps
    close(f) # input file

    # read reference xcorrs
    ref_f = jldopen(refname)

    # dictionary to keep track of convergence vectors for each station pair
    conv = Dict()
    # dictionary to keep track of stacked data
    stackData = Dict{String, Array{Float64,2}}()

    # iterate over timestamps
    for tstamp in tslist
        println("Processing $tstamp")
        # read xcorrs for each time step
        f_cur = jldopen(basefiname*".$tstamp.jld2")

        # iterate over station pairs
        for stnpair in keys(f_cur[tstamp])
            # load unstacked cross-correlations
            data = try f_cur["$tstamp/$stnpair"] catch; continue end

            # load reference cross-correlation for the current station pair
            ref = ref_f[stnpair].corr

            # iterate over windowed cross-correlations
            for i=1:size(data.corr)[2]
                if skipoverlap i=2*i-1 end # assumes cc_step is 1/2 cc_len

                # build input to convergence functions
                if haskey(stackData, stnpair) && metric=="rms"
                    # keep running sum of cross-correlations, stack if key exists, assign if not
                    stackData[stnpair] .+= data.corr[:, i]
                elseif haskey(stackData, stnpair) && metric=="snr"
                    # create matrix of xcorrs
                    stackData[stnpair] = hcat(stackData[stnpair], data.corr[:,i])
                else
                    # reshape into 2D array for snr input
                    stackData[stnpair] = reshape(data.corr[:, i], length(data.corr[:,i]), 1)
                end

                # compute convergence data point and append to vector
                if metric=="rms"
                    # compute RMS error of the running sum with the reference
                    if haskey(conv, stnpair)
                        push!(conv[stnpair], rms(normalize(ref[:,1]), normalize(stackData[stnpair][:,1])))
                    else
                        conv[stnpair] = [rms(normalize(ref[:,1]), normalize(stackData[stnpair][:,1]))]
                    end
                elseif metric=="snr"
                    # compute peak signal-to-noise ratio of increasingly long stacks
                    psnr = maximum(snr(stackData[stnpair], data.fs))
                    if haskey(conv, stnpair)
                        push!(conv[stnpair], psnr)
                    else
                        conv[stnpair] = [psnr]
                    end
                end
            end
        end
        close(f_cur) # xcorrs for each time step
    end
    close(ref_f) # reference xcorrs

    # write convergence vectors to disk
    f_conv = jldopen(foname, "a+")
    for key in keys(conv)
        f_conv[key] = conv[key]
    end
    close(f_conv) # convergence file
end
