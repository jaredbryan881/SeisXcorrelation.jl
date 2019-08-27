__precompile__()
#module SeisXcorrelation

using Distributed
@everywhere using SeisIO, SeisNoise, Dates, FFTW, JLD2, Sockets, ORCA
# using SeisIO, SeisNoise, Dates, FFTW, JLD2, Distributed, Sockets, ORCA

# read and write using jld2
@everywhere include("/n/home03/kokubo/.julia/dev/SeisXcorrelation/src/io.jl")
# generate, sort, and classify station pairs
@everywhere include("/n/home03/kokubo/.julia/dev/SeisXcorrelation/src/pairing.jl")
# miscellaneous utilities such as distance computation, windowing, and normalization
@everywhere include("/n/home03/kokubo/.julia/dev/SeisXcorrelation/src/utils.jl")
# splitting utilities for high-order cross-correlation coda extraction
@everywhere include("/n/home03/kokubo/.julia/dev/SeisXcorrelation/src/partition.jl")
# FFT wrapper for high-order cross-correlation
@everywhere include("/n/home03/kokubo/.julia/dev/SeisXcorrelation/src/fft.jl")
# cross-correlation wrapper for high-order cross-correlation
@everywhere include("/n/home03/kokubo/.julia/dev/SeisXcorrelation/src/correlate.jl")

export seisxcorrelation, seisxcorrelation_highorder


"""
    seisxcorrelation(tstamp::String, InputDict::Dict)

Compute cross-correlation save data in jld2 file with CorrData format.

# Arguments
- `data`::Dict    : !Reprecated! Dictionary containing input data in the form: tstamp/stn => SeisChannel
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `tserrorList`::Array{String,1}    : Array containing the name of failed cross-correlation timestamp/stationpairs
- `basefoname.tstamp.jld2`    : contains CorrData structure with a hierarchical structure (CC function, metadata)

## Update: Compute fft for the first time and store it into `FFTDict`. Xcorr from FFTDict afterwards, saving time
for loading SeisData and for compute_fft.
"""
function seisxcorrelation(tstamp::String, InputDict::Dict)
    # IO parameters
    finame     = InputDict["finame"]
    basefoname = InputDict["basefoname"]
    time_unit  = InputDict["timeunit"]
    # FFT parameters
    freqmin    = InputDict["freqmin"]
    freqmax    = InputDict["freqmax"]
    fs         = InputDict["fs"]
    cc_len     = InputDict["cc_len"]
    cc_step    = InputDict["cc_step"]
    to_whiten  = InputDict["to_whiten"]
    time_norm  = InputDict["time_norm"]
    maxdistance= InputDict["maxdistance"]
    # correlation parameters
    corrtype   = InputDict["corrtype"]   # "xcorr," "acorr," or "xchancorr" (or any combination thereof)
    corrmethod = InputDict["corrmethod"] # "cross-correlation", "deconv", or "coherence"
    maxtimelag = InputDict["maxtimelag"]
    # stacking parameters
    stack      = InputDict["allstack"]

    # dictionary to cache FFTs
    FFTDict = Dict{String, FFTData}()
    # list of stations that failed for this timestep
    tserrorList = []

    # load input file
    inFile = jldopen(finame, "r")

    # split dataset names (keys of data) by "/" to get station list
    # assume form "$tstamp/$station"
    #dsk = collect(keys(data))
    #stlist = sort([string(split.(dsk[i], "/")[2]) for i=1:length(dsk)])
    stlist = collect(keys(inFile["$tstamp"]))

    # create output file for each time stamp, fill relevant info
    outFile = jldopen("$basefoname.$tstamp.jld2", "a+")
    outFile["info/stationlist"] = stlist
    outFile["info/timeunit"] = time_unit

    println("$tstamp: Computing cross-correlations")
    stniter = 0 # counter to prevent computing duplicate xcorrs with reversed order
    # iterate over station list

    for stn1 in stlist
        stniter+=1

        # don't attempt FFT if this failed already
        if stn1 in tserrorList continue end

        FFT1 = try
            # try to read FFT from cached FFTs
            if stn1 in keys(FFTDict)
                FFTDict[stn1]
            else

                # read station SeisChannels into SeisData before FFT
                S1 = try
                    #SeisData(data["$tstamp/$stn1"])
                    SeisData(inFile["$tstamp/$stn1"])
                catch;
                    push!(tserrorList, "$tstamp/$stn1")
                    filter!(a->a≠stn1, stlist)
                    println("$tstamp: $stn1 encountered an error on read. Skipping.")
                    continue
                end

                if size(S1)[1] > 1 @warn "SeisData contains multiple channels. Operating only on the first." end
                delete!(S1[1].misc, "kurtosis")
                delete!(S1[1].misc, "eqtimewindow")

                # do not attempt fft if data was not available
                try
                    if S1[1].misc["dlerror"] == 1
                        push!(tserrorList, "$tstamp/$stn1")
                        filter!(a->a≠stn1, stlist)
                        println("$tstamp: $stn1 encountered an error: dlerror==1. Skipping.")
                        continue
                    end
                catch;
                    # assume key "dlerror" does not exist (not downloaded via SeisDownload)
                end

                # make sure the data is the proper length to avoid dimension mismatch
                npts1 = Int(time_unit * S1[1].fs)
                if (length(S1[1].x) > npts1) S1[1].x=S1[1].x[1:npts1]; S1[1].t[end,1]=npts1 end

                tt1temp = @elapsed FFT1 = compute_fft(S1, freqmin, freqmax, fs, Int(cc_step), Int(cc_len), to_whiten=to_whiten, time_norm=time_norm)
                #println("fft1: $tt1temp ")
                FFTDict[stn1] = FFT1
                FFT1
            end
        catch y
            push!(tserrorList, "$tstamp/$stn1")
            filter!(a->a≠stn1, stlist)
            println("$tstamp: $stn1 encountered an error on FFT1. Skipping.")
            continue
        end

        # Evaluate Memory use
        # if stniter == 1
        #     memory = trunc(Int64, sizeof(FFT1) * length(stlist) / 1e3) #[GB]
        #     println("Maximum memory use can be: $memory [KB]")
        # end

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

                # try to read FFT from cached FFTs
                FFT2 = try
                    if stn2 in keys(FFTDict)
                        FFTDict[stn2]
                    else
                        # read station SeisChannels into SeisData before FFT
                        S2 = try
                            #SeisData(data["$tstamp/$stn2"])
                            SeisData(inFile["$tstamp/$stn2"])
                        catch;
                            push!(tserrorList, "$tstamp/$stn2")
                            filter!(a->a≠stn2, stlist)
                            println("$tstamp: $stn2 encountered an error on read. Skipping.")
                            continue
                        end

                        if size(S2)[1] > 1 @warn "SeisData contains multiple channels. Operating only on the first." end
                        delete!(S2[1].misc, "kurtosis")
                        delete!(S2[1].misc, "eqtimewindow")

                        # do not attempt fft if data was not available
                        try
                            if S2[1].misc["dlerror"] == 1
                                push!(tserrorList, "$tstamp/$stn2")
                                filter!(a->a≠stn2, stlist)
                                println("$tstamp: $stn2 encountered an error: dlerror==1. Skipping.")
                                continue
                            end
                        catch;
                            # assume that the key dlerror does not exist (data not downloaded via SeisDownload)
                        end

                        # make sure the data is the proper length to avoid dimension mismatch
                        npts2 = Int(time_unit * S2[1].fs)
                        if (length(S2[1].x) > npts2) S2[1].x=S2[1].x[1:npts2]; S2[1].t[end,1]=npts2 end

                        tt2temp = @elapsed FFT2 = compute_fft(S2, freqmin, freqmax, fs, Int(cc_step), Int(cc_len), to_whiten=to_whiten, time_norm=time_norm)
                        #print("fft2: $tt2temp ")
                        FFTDict[stn2] = FFT2
                        FFT2
                    end
                catch y
                    push!(tserrorList, "$tstamp/$stn2")
                    filter!(a->a≠stn2, stlist)
                    println("$tstamp: $stn2 encountered an error on FFT2. Skipping.")
                    continue
                end

            else
                #println("Skipping cross-correlation of $stn1 and $stn2.")
                continue
            end

            # do not attempt fft if distance of station pair is too large (dist > maxdistance)
            #print(dist(S1[1].loc, S2[1].loc)./1e3)
            #println(" [km]")
            if maxdistance == false
                # compute cc
            elseif dist(FFT1.loc, FFT2.loc) > maxdistance
                #println("$tstamp: distance $stn1 - $stn2 is greater than maxdistance. Skipping.")
                continue;
            end

            # compute correlation using SeisNoise.jl -- returns type CorrData
            tt3temp = @elapsed xcorr = compute_cc(FFT1, FFT2, maxtimelag, corr_type=corrmethod)
            #print("xcorr: $tt3temp ")
            varname = "$tstamp/$stn1.$stn2"
            try
                # compute distance between stations
                xcorr.misc["dist"] = dist(FFT1.loc, FFT2.loc)
                # save location of each station
                xcorr.misc["location"] = Dict(stn1=>FFT1.loc, stn2=>FFT2.loc)

                # stack over DL_time_unit
                if stack==true stack!(xcorr, allstack=true) end

                outFile[varname] = xcorr
            catch
                println("$stn1 and $stn2 have no overlap at $tstamp.")
                push!(tserrorList, "$tstamp/$stn1")
            end

        end

        # release memory held by the FFT and time series for this station
        delete!(FFTDict, stn1)
    end
    outFile["info/errors"] = tserrorList
    close(inFile)
    close(outFile)
end

"""

    seisxcorrelation_highorder(tstamp::String, corrstationlist::Array{String,2}, data::JLD2.JLDFile, InputDict::Dict)

Compute C2 or C3 cross-correlation and save data in jld2 file with CorrData format.

# Arguments
- `tstamp::String,`    : Time stamp to read from JLD2 file.
- `corrstationlist::Array{String,1},`    : List of cross-correlation names, e.g. BP.CCRB..BP1.BP.EADB..BP1
- `data::JLD2.JLDFile`    : JLD2 file object to read data from.
- `InputDict::Dict`    : Dictionary containing IO, FFT, xcorr, and stacking parameters

# Output
- `basefoname.tstamp.jld2`    : contains SeisData structure with a hierarchical structure (CC function, metadata)

"""
function seisxcorrelation_highorder(data::Dict, tstamp::String, corrstationlist::Array{String,2}, allowable_pairs::Array{String,1}, allowable_vs::Array{String,1}, InputDict::Dict)
    # IO parameters
    basefoname = InputDict["basefoname"]
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
    maxtimelag = InputDict["maxtimelag"]
    # stacking parameters
    stack      = InputDict["allstack"]

    pairiter = 1 # counter to prevent computing duplicate xcorrs with reversed order

    # dictionary to cache FFTs
    FFTDict = Dict{String, Array{FFTData,1}}()

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

    outFile = jldopen("$basefoname.$tstamp.jld2", "a+")
    outFile["info/corrstationlist"] = xcorrlist

    println("$tstamp: Computing cross-correlations")
    for p1=1:length(xcorrlist[1,:])
        # first station pair name and its reverse
        pair1 = xcorrlist[:, p1]
        p1name = pair1[1]*"."*pair1[2]
        p1namerev = pair1[2]*"."*pair1[1]

        for p2=pairiter:length(xcorrlist[1,:])
            # second station pair name
            pair2 = xcorrlist[:, p2]
            p2name = pair2[1]*"."*pair2[2]

            # read data only if we have 3 unique stations
            virt_src = intersect(pair1, pair2)

            if length(virt_src) != 1 continue end

            # set stations used as receivers for c3
            stn1 = setdiff(pair1, virt_src)[1]
            stn2 = setdiff(pair2, virt_src)[1]

            # keep string instead of array for ease of use
            virt_src = virt_src[1]

            # read data only if we have 3 unique station pairs
            if get_corrtype([stn1, stn2]) != "xcorr" continue end
            # read data only if the virtual source is approved
            if virt_src ∉ allowable_vs continue end

            # read data only if the station pair is approved
            if ("$stn1.$stn2" ∉ allowable_pairs) && ("$stn2.$stn1" ∉ allowable_pairs) continue end

            # load FFT1 from dict if it exists
            if "$(virt_src)/$stn1" in keys(FFTDict)
                FFT1_neg, FFT1_pos = FFTDict["$(virt_src)/$stn1"]
            else
                # load data and (optionally) reverse it to get virt_src-stn config
                # copy struct to prevent in-place reversal from messing up future xcorrs
                xcorr1 = copy(data["$tstamp/$p1name"])

                if size(xcorr1.corr, 2)>1 stack!(xcorr1) end

                # reverse xcorr if we loaded stn1-virt_src to get virt_src-stn1
                if virt_src==pair1[2]
                    xcorr1.corr=reverse(xcorr1.corr, dims=1)
                    vs_loc = xcorr1.misc["location"][pair1[2]]
                    stn1_loc = xcorr1.misc["location"][pair1[1]]
                else
                    vs_loc = xcorr1.misc["location"][pair1[1]]
                    stn1_loc = xcorr1.misc["location"][pair1[2]]
                end

                # partition xcorr into positive and negative windows
                xcorr1.corr = partition(xcorr1.corr, start_lag, win_len)

                # compute FFT
                FFT1_neg, FFT1_pos = compute_fft_c3(xcorr1, freqmin, freqmax, fs, cc_step, cc_len)

                FFTDict["$virt_src/$stn1"] = [FFT1_neg, FFT1_pos]
            end

            # load FFT2 from dict if it exists
            if "$virt_src/$stn2" in keys(FFTDict)
                FFT2_neg, FFT2_pos = FFTDict["$virt_src/$stn2"]
            else
                # load data and (optionally) reverse it to get virt_src-stn config
                # copy struct to prevent in-place reversal from messing up future xcorrs
                xcorr2 = copy(data["$tstamp/$p2name"])

                if size(xcorr2.corr, 2)>1 stack!(xcorr2) end

                # reverse xcorr if we loaded stn1-virt_src to get virt_src-stn1
                if virt_src==pair2[2]
                    xcorr2.corr=reverse(xcorr2.corr, dims=1)
                    stn2_loc = xcorr2.misc["location"][pair2[1]]
                else
                    stn2_loc = xcorr2.misc["location"][pair2[2]]
                end

                # partition xcorr into positive and negative windows
                xcorr2.corr = partition(xcorr2.corr, start_lag, win_len)

                # compute FFT
                FFT2_neg, FFT2_pos = compute_fft_c3(xcorr2, freqmin, freqmax, fs, cc_step, cc_len)

                FFTDict["$virt_src/$stn2"] = [FFT2_neg, FFT2_pos]
            end

            #compute xcorr - 1st half of dim=2 is neg lag, 2nd half is pos lag
            xcorr_c3_neg = compute_cc_c3(FFT1_neg, FFT2_neg, maxtimelag)
            xcorr_c3_neg.name = "$stn1.$stn2.$(virt_src)_neg"
            xcorr_c3_pos = compute_cc_c3(FFT1_pos, FFT2_pos, maxtimelag)
            xcorr_c3_pos.name = "$stn1.$stn2.$(virt_src)_pos"

            # stacking
            if stack
                stack!(xcorr_c3_neg, allstack=true)
                stack!(xcorr_c3_pos, allstack=true)
            end

            # save cross-correlation
            varname = "$tstamp/$stn1.$stn2/$virt_src"
            try
                outFile[varname] = [xcorr_c3_neg, xcorr_c3_pos]
            catch
                println("$varname encountered an error on saving to JLD2.")
                push!(tserrorList, varname)
                continue
            end
        end

        # release memory held by the FFT for this station pair in FFTDict
        delete!(FFTDict, p1name)
        delete!(FFTDict, p1namerev)
        # release memory held by the input data
        delete!(data, "$tstamp/$p1name")
        delete!(data, "$tstamp/$p1namerev")

        pairiter += 1
    end
    outFile["erorrs"] = tserrorList
    close(outFile)
end

#end # module
