__precompile__()
#module SeisXcorrelation

using SeisIO, SeisDownload, Noise, Printf, Dates, FFTW, JLD2, Distributed

include("pairing.jl")
include("fft.jl")
include("utils.jl")
#export seisxcorrelation


"""

    seisxcorrelation(maxtimelag::Real, finame::String, foname::String; IsAllComponent::Bool=false)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `finame::String,`    : Input file name e.g. network = "BPnetwork"
- `maxtimelag::Real`    : Maximum lag time e.g. maxtimelag = 100 [s]
- `corrtype::AbstractArray`    : Array of strings containing types of correlations to compute, e.g. ["xcorr", "acorr"]
- `corrorder`::Int    : Order of cross correlation, e.g. 3 for C3 (high order cross correlation)
- `IsAllComponent::Bool`   : If true, compute 3 X 3 components cross-correlation

# Output
- `foname.jld2`    : contains SeisData structure with a hierarchical structure (CC function, metadata)

"""
function seisxcorrelation(tstamp::String, InputDict::Dict)
    # IO parameters
    finame     = InputDict["finame"]
    foname     = InputDict["foname"]
    # correlation parameters
    corrtype   = InputDict["corrtype"]
    corrorder  = InputDict["corrorder"]
    maxtimelag = InputDict["maxtimelag"]
    # FFT parameters
    freqmin    = InputDict["freqmin"]
    freqmax    = InputDict["freqmax"]
    fs         = InputDict["fs"]
    cc_len     = InputDict["cc_len"]
    cc_step    = InputDict["cc_step"]
    # stacking parameters
    stack      = InputDict["allstack"]

    stniter = 1 # counter to enforce computing only unique cross-correlations

    # dictionary to contain errors
    errorDict = Dict{String, Array{String,1}}("DimensionMismatch"=>[""], "DataUnavailable"=>[""])
    # dictionary to cache FFTs
    FFTDict = Dict{String, FFTData}()
    # list of stations that failed for this timestep
    tserrorList = []

    # iterate over station list
    for stn1 in stlist
        if "$stn1" in tserrorList continue end # don't attempt FFT if this failed already
        # read station SeisChannels into SeisData before FFT
        S1 = SeisData(data["$tstamp/$stn1"])

        # do not attempt fft if data was not available
        if S1[1].misc["dlerror"] == 1
            push!(errorDict["DataUnavailable"], "$stn1")
            push!(tserrorList, "$tstamp/$stn1")
            println("$tstamp/$stn1 unavailable.")
            continue
        end

        # round start times to nearest millisecond to avoid split start times bug
        S1[1].t[1,2] = round(S1[1].t[1,2], sigdigits=13)

        # check correlation order and compute the appropriate FFT
        if corrorder == 1
            FFT1 = try
                # try to read FFT from cached FFTs
                if "$stn1" in keys(FFTDict)
                    FFTDict["$stn1"]
                else
                    FFT1 = compute_fft(S1, freqmin, freqmax, fs, cc_step, cc_len)
                    FFTDict["$stn1"] = FFT1
                    FFT1
                end
            catch y
                if isa(y, DimensionMismatch)
                    # TODO handle dimension mismatch in a way that doesnt require throwing away the full day
                    push!(errorDict["DimensionMismatch"], "$stn1")
                    push!(tserrorList, "$tstamp/$stn1")
                    println("$tstamp: $stn1 has a dimension mismatch")
                    continue
                end
            end
        else
            println("Higher order correlation methods are not yet implemented\n")
            exit()
        end

        # iterate over station list again
        for stn2 in stlist[stniter:end]
            if "$stn2" in tserrorList continue end # don't attempt FFT if this failed already
            println("$tstamp: corrrelating $stn1 with $stn2")

            # see if this is an auto-, cross-, or xchan-correlation
            ct = get_corrtype([stn1, stn2])

            # autocorrelation
            if (ct=="acorr") && ("acorr" in corrtype)
                # set the stn2 FFT to the already computed FFT for stn1
                FFT2 = FFT1

            # cross-channel correlation
            elseif (ct=="xchancorr") && ("xchancorr" in corrtype)
                # read station SeisChannels into SeisData before FFT
                S2 = SeisData(data["$tstamp/$stn2"])

                # do not attempt fft if data was not available
                if S2[1].misc["dlerror"] == 1
                    push!(errorDict["DataUnavailable"], "$stn2")
                    push!(tserrorList, "$tstamp/$stn1")
                    println("$tstamp/$stn2 unavailable.")
                    continue
                end

                # round start times to nearest millisecond to avoid split start times bug
                S2[1].t[1,2] = round(S2[1].t[1,2], sigdigits=13)
                # check correlation order and compute the appropriate FFT using Noise.jl
                if corrorder == 1
                    # try to read FFT from cached FFTs
                    FFT2 = try
                        if "$stn2" in keys(FFTDict)
                            FFTDict["$stn2"]
                        else
                            FFT2 = compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len)
                            FFTDict["$stn2"] = FFT2
                            FFT2
                        end
                    catch y
                        if isa(y, DimensionMismatch)
                            push!(errorDict["DimensionMismatch"], "$stn2")
                            push!(tserrorList, "$tstamp/$stn2")
                            println("$tstamp: $stn2 has a dimension mismatch")
                            continue
                        end
                    end
                else
                    println("Higher order correlation methods are not yet implemented\n")
                    exit()
                end

            # cross-correlation
            elseif (ct=="xcorr") && ("xcorr" in corrtype)
                # read station SeisChannels into SeisData before FFT
                S2 = SeisData(data["$tstamp/$stn2"])

                # do not attempt fft if data was not available
                if S2[1].misc["dlerror"] == 1
                    push!(errorDict["DataUnavailable"], "$stn2")
                    push!(tserrorList, "$tstamp/$stn1")
                    println("$tstamp/$stn2 unavailable.")
                    continue
                end

                # round start times to nearest millisecond to avoid split start times bug
                S2[1].t[1,2] = round(S2[1].t[1,2], sigdigits=13)

                # check correlation order and compute the appropriate FFT using Noise.jl
                if corrorder == 1
                    # try to read FFT from cached FFTs
                    FFT2 = try
                        if "$stn2" in keys(FFTDict)
                            FFTDict["$stn2"]
                        else
                            FFT2 = compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len)
                            FFTDict["$stn2"] = FFT2
                            FFT2
                        end
                    catch y
                        if isa(y, DimensionMismatch)
                            push!(errorDict["DimensionMismatch"], "$stn2")
                            push!(tserrorList, "$tstamp/$stn2")
                            println("$tstamp: $stn2 has a dimension mismatch")
                            continue
                        end
                    end
                else corrorder == 2
                    println("Higher order correlation methods are not yet implemented\n")
                    exit()
                end

            else
                println("Skipping cross-correlation of $stn1 and $stn2.")
                continue
            end

            # compute correlation using Noise.jl -- returns type CorrData
            xcorr = compute_cc(FFT1, FFT2, maxtimelag)

            if stack==true
                # stack hourly cross-correlations for a single daily average
                stack!(xcorr, allstack=true)
            end

            # save data after each cross correlation
            varname = "$tstamp/$stn1.$stn2"
            try
                save_CorrData2JLD2(foname, varname, xcorr)
            catch
                println("$stn1 and $stn2 have no overlap at $tstamp.")
                push!(tserrorList, "$tstamp/$stn1")
            end
        end
        stniter += 1
    end
    # save dict of dim mismatch and data availability errors to JLD2
    save_Dict2JLD2(foname, "$tstamp/ccerrors", errorDict)
    # save list of all cross correlation names that failed
    return tserrorList
end
