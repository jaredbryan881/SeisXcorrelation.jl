using SeisIO, JLD2, Noise, PlotlyJS
include("utils.jl")

"""

    compute_reference_xcorr(finame::String, foname::String)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `finame::String,`    : Input file name e.g. "./inputData/BPnetwork.jld2"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"
- `xcorrSize::Int,`    : Number of points in a cross-correlation window

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair

"""
function compute_reference_xcorr(finame::String, foname::String, xcorrSize::Int)
    f = jldopen(finame)

    corrstationlist = f["info/corrstationlist"] # names of station pairs
    tslist = f["info/timestamplist"] # time stamps
    reference_xcorrnames = []

    # iterate over correlation type
    for ct in keys(corrstationlist)
        # iterate over station pairs
        for stnum=1:size(corrstationlist[ct])[2]
            # get name of station pair
            stn1 = corrstationlist[ct][1, stnum]
            stn2 = corrstationlist[ct][2, stnum]
            varname = "$stn1" .* "." .* "$stn2"
            println(varname)

            # initialize array to hold stacked data
            stackData = zeros(xcorrSize)
            # num of summed xcorrs
            numxcorr = 0
            # compute full stack for this station pair
            for tstamp in tslist
                # do not attempt stacking if there was an error in cross-correlation
                if "$tstamp/$stn1" in f["info/tserrors"] || "$tstamp/$stn2" in f["info/tserrors"]
                    continue
                end

                # read data and stack
                data = try f["$tstamp/$varname"].corr catch; continue end
                stackData += sum(data, dims=2)
                numxcorr += size(data)[2]
            end
            # compute average
            stackData /= numxcorr

            # save reference xcorr
            save_Array2JLD2(foname, "$varname", stackData)
            push!(reference_xcorrnames, varname)
        end
    end
    close(f)
    save_Array2JLD2(foname, "info/reference_xcorrnames", reference_xcorrnames)
end

"""

    convergence(finame::String, foname::String)

Compute the RMS error between progressively longer stacks of cross-correlations and a reference cross-correlation

# Arguments
- `finame::String,`    : Input file name for unstacked cross-correlations e.g. "./inputData/BPnetwork.jld2"
- `finame::refname,`    : Input file name for reference cross-corrleations e.g. "./inputData/BPnetwork.jld2"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"
- `ntwindows::Int`    : Number of cross-correlations in each CorrData object
- `skipoverlap::Bool`    : Whether or not to skip overlap in cross-correlation windows. If cross-correlations were computed
                           with no overlap, set this to false to ensure no time is skipped.

# Output
- `foname.jld2`    : contains arrays of RMS errors for progressively longer stacks of cross-correlations

"""
function convergence(finame::String, refname::String, foname::String, nxcorr::Int, skipoverlap::Bool=false)
    f = jldopen(finame)
    ref_f = jldopen(refname)

    tslist = f["info/timestamplist"] # time stamps
    stname = ref_f["info/reference_xcorrnames"] # station pairs with computed reference cross-correlations

    # iterate over station pairs
    for stpair in stname
        # reference xcorr
        ref = ref_f[stpair]
        # initialize running stack with size of the reference
        stackData = zeros(size(ref))
        # initialize rmse array
        rms_conv = zeros(length(tslist)*nxcorr)

        iter=1
        # iterate over timestamps
        for tstamp in tslist
            data = try f["$tstamp/$stpair"].corr catch; continue end
            # iterate over windowed cross-correlations
            for i=1:nxcorr
                if skipoverlap
                    # assume overlap is half of window length and skip overlapping windows
                    i=2*i - 1
                end
                # keep running sum of cross-correlations
                stackData += data[:, i]

                # compute RMS error of the average of the running sum with the reference
                rms_conv[iter] = rms(ref[:,1], stackData[:,1] ./ iter)
                # keep track of the number of xcorrs in the running sum
                iter += 1
            end
        end
        # write RMS errors to foname
        save_Array2JLD2(foname, stpair, rms_conv)
    end
end

#compute_reference_xcorr("../EXAMPLE/xcorr_BP/outputData/BPnetworkxcorr_weq.jld2", "reference_xcorr.jld2", 4001)
convergence("../EXAMPLE/xcorr_BP/outputData/BPnetworkxcorr_weq.jld2", "reference_xcorr.jld2", "rms_weq_wol.jld2", 46)
