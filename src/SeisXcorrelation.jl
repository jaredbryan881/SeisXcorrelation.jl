__precompile__()
#module SeisXcorrelation

using SeisIO, SeisNoise, Dates, FFTW, JLD2, Distributed, PlotlyJS, Sockets, ORCA

include("pairing.jl")
include("fft.jl")
include("utils.jl")
include("io.jl")
include("partition.jl")
include("correlate.jl")

export seisxcorrelation, seisxcorrelation_highorder


"""

    seisxcorrelation(tstamp::String, stationlist::Array{String,1}, inputData::JLD2.JLDFile, InputDict::Dict)

Compute cross-correlation save data in jld2 file with CorrData format.

# Arguments
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `stationlist::Array{String,1}`    : List of stations to be cross-correlation
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `tserrorList`::Array{String,1}    : Array containing the name of failed cross-correlation timestamp/stationpairs
- `basefoname.tstamp.jld2`    : contains SeisData structure with a hierarchical structure (CC function, metadata)

"""
function seisxcorrelation(data::Dict, tstamp::String, stationlist::Array{String,1}, InputDict::Dict)
    # IO parameters
    basefoname = InputDict["basefoname"]
    time_unit  = InputDict["timeunit"]
    # FFT parameters
    freqmin    = InputDict["freqmin"]
    freqmax    = InputDict["freqmax"]
    fs         = InputDict["fs"]
    cc_len     = InputDict["cc_len"]
    cc_step    = InputDict["cc_step"]
    # correlation parameters
    corrtype   = InputDict["corrtype"]
    corrorder  = InputDict["corrorder"]
    maxtimelag = InputDict["maxtimelag"]
    # stacking parameters
    stack      = InputDict["allstack"]

    stniter = 1 # counter to prevent computing duplicate xcorrs with reversed order

    # dictionary to cache FFTs
    FFTDict = Dict{String, FFTData}()

    # list of stations that failed for this timestep
    tserrorList = []

    # copy station list to avoid problems when removing stations
    stlist = deepcopy(stationlist)

    # create output file for each time stamp, fill relevant info
    outFile = jldopen("$basefoname.$tstamp.jld2", "a+")
    outFile["info/stationlist"] = stlist
    outFile["info/timeunit"] = time_unit

    println("$tstamp: Computing cross-correlations")
    # iterate over station list
    for stn1 in stlist
        # don't attempt FFT if this failed already
        if stn1 in tserrorList continue end

        # read station SeisChannels into SeisData before FFT
        S1 = SeisData(data["$tstamp/$stn1"])

        # do not attempt fft if data was not available
        if S1[1].misc["dlerror"] == 1
            push!(tserrorList, "$tstamp/$stn1")
            filter!(a->a≠stn1, stlist)
            println("$tstamp: $stn1 encountered an error. Skipping.")
            continue
        end

        # make sure the data is the proper length to avoid dimension mismatch
        numpts1 = time_unit * S1[1].fs
        if (length(S1[1].x) > numpts1) S1[1].x=S1[1].x[1:npts1]; S1[1].t[2,1]=npts1 end

        FFT1 = try
            # try to read FFT from cached FFTs
            if stn1 in keys(FFTDict)
                FFTDict[stn1]
            else
                FFT1 = compute_fft(S1, freqmin, freqmax, fs, cc_step, cc_len)
                FFTDict[stn1] = FFT1
                FFT1
            end
        catch y
            push!(tserrorList, "$tstamp/$stn1")
            filter!(a->a≠stn1, stlist)
            println("$tstamp: $stn1 encountered an error. Skipping.")
            continue
        end

        # iterate over station list again
        for stn2 in stlist[stniter:end]
            if stn2 in tserrorList continue end # don't attempt FFT if this failed already

            # see if this is an auto-, cross-, or xchan-correlation
            ct = get_corrtype([stn1, stn2])

            # autocorrelation
            if (ct=="acorr") && ("acorr" in corrtype)
                # set the stn2 FFT to the already computed FFT for stn1
                FFT2 = FFT1

            # cross- or cross-channel correlation
            elseif ct in corrtype
                # read station SeisChannels into SeisData before FFT
                S2 = SeisData(data["$tstamp/$stn2"])

                # do not attempt fft if data was not available
                if S2[1].misc["dlerror"] == 1
                    push!(tserrorList, "$tstamp/$stn2")
                    filter!(a->a≠stn2, stlist)
                    println("$tstamp: $stn2 encountered an error. Skipping.")
                    continue
                end

                # make sure the data is the proper length to avoid dimension mismatch
                numpts2 = time_unit * S2[1].fs
                if (length(S2[1].x) > numpts2) S2[1].x=S2[1].x[1:npts2]; S2[1].t[2,1]=npts2 end

                # try to read FFT from cached FFTs
                FFT2 = try
                    if stn2 in keys(FFTDict)
                        FFTDict[stn2]
                    else
                        FFT2 = compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len)
                        FFTDict[stn2] = FFT2
                        FFT2
                    end
                catch y
                    push!(tserrorList, "$tstamp/$stn2")
                    filter!(a->a≠stn2, stlist)
                    println("$tstamp: $stn2 encountered an error. Skipping.")
                    continue
                end

            else
                println("Skipping cross-correlation of $stn1 and $stn2.")
                continue
            end

            # compute correlation using Noise.jl -- returns type CorrData
            xcorr = compute_cc(FFT1, FFT2, maxtimelag)
            # compute distance between stations
            xcorr.misc["dist"] = dist(FFT1.loc, FFT2.loc)
            # TODO: make sure both station locations are maintained in CorrData
            # stack over DL_time_unit
            if stack==true stack!(xcorr, allstack=true) end

            varname = "$tstamp/$stn1.$stn2"
            try
                outFile[varname] = xcorr
            catch
                println("$stn1 and $stn2 have no overlap at $tstamp.")
                push!(tserrorList, "$tstamp/$stn1")
            end
        end
        stniter += 1

        # release memory held by the FFT and time series for this station
        delete!(FFTDict, stn1)
        delete!(data, "$tstamp/$stn1")
    end

    close(outFile)
    return tserrorList
end

"""

    seisxcorrelation_highorder(tstamp::String, corrstationlist::Array{String,2}, inputData::JLD2.JLDFile, InputDict::Dict)

Compute C2 or C3 cross-correlation and save data in jld2 file with CorrData format.

# Arguments
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `corrstationlist::Array{String,1},`    : List of cross-correlation names, e.g. BP.CCRB..BP1.BP.EADB..BP1
- `inputData::JLD2.JLDFile`    : JLD2 file object to read data from.
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `foname.jld2`    : contains SeisData structure with a hierarchical structure (CC function, metadata)

"""
function seisxcorrelation_highorder(tstamp::String, corrstationlist::Array{String,2}, allowable_pairs::Array{String,1}, allowable_vs::Array{String,1}, inputData::JLD2.JLDFile, InputDict::Dict)
    # IO parameters
    foname     = InputDict["foname"]
    # Partitioning parameters
    start_lag  = InputDict["start_lag"]
    win_len    = InputDict["window_len"]
    # FFT parameters
    freqmin    = InputDict["freqmin"]
    freqmax    = InputDict["freqmax"]
    fs         = InputDict["fs"]
    cc_len     = InputDict["cc_len"]
    cc_step    = InputDict["cc_step"]
    # correlation parameters
    corrtype   = InputDict["corrtype"]
    corrorder  = InputDict["corrorder"]
    maxtimelag = InputDict["maxtimelag"]
    # stacking parameters
    stack      = InputDict["allstack"]

    pairiter = 1 # counter to prevent computing duplicate xcorrs with reversed order

    # dictionary to cache FFTs
    FFTDict = Dict{String, FFTData}()

    # list of stations that failed for this timestep
    tserrorList = []

    xcorrlist = deepcopy(corrstationlist)

    # convert samples given in [s] to the number of samples
    if typeof(start_lag) == Float64
        start_lag = convert(Int64, start_lag*fs)
    end
    if typeof(win_len) == Float64
        win_len  = convert(Int64, win_len*fs)
    end

    for p1=1:length(xcorrlist[1,:])
        pair1 = xcorrlist[:, p1]
        p1name = pair1[1]*"."*pair1[2]
        for p2=pairiter:length(xcorrlist[1,:])
            pair2 = xcorrlist[:, p2]
            p2name = pair2[1]*"."*pair2[2]

            # read data only if we have 3 unique stations
            virt_src = intersect(pair1, pair2)
            if length(virt_src) != 1 continue end

            # set stations used as receivers for c3
            stn1 = setdiff(pair1, virt_src)[1]
            stn2 = setdiff(pair2, virt_src)[1]

            # read data only if we have 3 unique station pairs
            if get_corrtype([stn1, stn2]) != "xcorr" continue end
            # read data only if the virtual source is approved
            if virt_src ∉ allowable_vs continue end
            # read data only if the station pair is approved
            if ("$stn1.$stn2" ∉ allowable_pairs) || ("$stn2.$stn1" ∉ allowable_pairs) continue end

            # load FFT1 from dict if it exists
            if "$(virt_src[1])/$stn1" in keys(FFTDict)
                FFT1 = FFTDict["$(virt_src[1])/$stn1"]
            else
                # load data and (optionally) reverse it to get virt_src-stn config
                # copy struct to prevent in-place reversal from messing up future xcorrs
                xcorr1 = copy(inputData["$tstamp/$p1name"])

                # reverse xcorr if we loaded stn1-virt_src to get virt_src-stn1
                if virt_src==pair1[2] xcorr1.corr=reverse(xcorr1.corr) end
                # partition xcorr into positive and negative windows
                xcorr1.corr = partition(xcorr1.corr, start_lag, win_len)

                # compute FFT
                FFT1 = compute_fft_c3(xcorr1, freqmin, freqmax, fs, cc_step, cc_len)

                FFTDict["$(virt_src[1])/$stn1"] = FFT1
            end

            # load FFT2 from dict if it exists
            if "$(virt_src[1])/$stn2" in keys(FFTDict)
                FFT2 = FFTDict["$(virt_src[1])/$stn2"]
            else
                # load data and (optionally) reverse it to get virt_src-stn config
                # copy struct to prevent in-place reversal from messing up future xcorrs
                xcorr2 = copy(inputData["$tstamp/$p2name"])

                # reverse xcorr if we loaded stn1-virt_src to get virt_src-stn1
                if virt_src==pair2[2] xcorr2.corr=reverse(xcorr2.corr) end
                # partition xcorr into positive and negative windows
                xcorr2.corr = partition(xcorr2.corr, start_lag, win_len)

                # compute FFT
                FFT2 = compute_fft_c3(xcorr2, freqmin, freqmax, fs, cc_step, cc_len)

                FFTDict["$virt_src/$stn2"] = FFT2
            end

            #compute xcorr - 1st half of dim=2 is neg lag, 2nd half is pos lag
            xcorr_c3 = compute_cc_c3(FFT1, FFT2, maxtimelag)
            xcorr_c3.name = "$stn1.$stn2.$(virt_src[1])"

            # save cross-correlation
            varname = "$tstamp/$stn1.$stn2/$(virt_src[1])"
            try
                save_CorrData2JLD2(foname, varname, xcorr)
            catch
                println("$varname encountered an error on saving to JLD2.")
                push!(tserrorList, varname)
                continue
            end
        end

        # release memory held by the FFT for this station pair in FFTDict
        delete!(FFTDict, p1name)
        delete!(FFTDict, pair1[2]*"."*pair1[1])

        pairiter += 1
    end

    return tserrors
end

#end # module
