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
function seisxcorrelation(tstamp::String, finame::String, foname::String, corrtype::Array{String,1}, corrorder::Int, maxtimelag::Real, freqmin::Real, freqmax::Real, fs::Real, cc_len::Int, cc_step::Int, IsAllComponents::Bool=false)
    # iterate over station list
    for stn1 in stlist
        # read station SeisChannels into SeisData before FFT
        S1 = SeisData(data["$tstamp/$stn1"])
        # round start times to nearest millisecond to avoid split start times bug
        S1[1].t[1,2] = round(S1[1].t[1,2], sigdigits=13)

        # check correlation order and compute the appropriate FFT
        if corrorder == 1
            FFT1 = compute_fft(S1, freqmin, freqmax, fs, cc_step, cc_len)
        else
            println("Higher order correlation methods are not yet implemented\n")
            exit()
        end

        # iterate over station list again
        for stn2 in stlist
            @time begin
            println("$tstamp: corrrelating $stn1 with $stn2")
            # see if this is an auto-, cross-, or xchan-correlation
            ct = get_corrtype([stn1, stn2])

            # check if this is an autocorrelation
            if (ct=="acorr") && ("acorr" in corrtype)
                # set the stn2 FFT to the already computed FFT for stn1
                FFT2 = FFT1

                # check if this is a cross-channel correlation
            elseif (ct=="xchancorr") && ("xchancorr" in corrtype)
                # read station SeisChannels into SeisData before FFT
                S2 = SeisData(data["$tstamp/$stn2"])
                # round start times to nearest millisecond to avoid split start times bug
                S2[1].t[1,2] = round(S2[1].t[1,2], sigdigits=13)
                # check correlation order and compute the appropriate FFT using Noise.jl
                if corrorder == 1
                    FFT2 = compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len)
                else
                    println("Higher order correlation methods are not yet implemented\n")
                    exit()
                end

                # check if this is a cross-correlation
            elseif (ct=="xcorr") && ("xcorr" in corrtype)
                # read station SeisChannels into SeisData before FFT
                S2 = SeisData(data["$tstamp/$stn2"])
                # round start times to nearest millisecond to avoid split start times bug
                S2[1].t[1,2] = round(S2[1].t[1,2], sigdigits=13)

                # check correlation order and compute the appropriate FFT using Noise.jl
                if corrorder == 1
                    FFT2 = compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len)
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

            # save data after each cross correlation
            varname = "$tstamp/$stn1.$stn2"
            try save_CorrData2JLD2(foname, varname, xcorr) catch; println("$stn1 and $stn2 have no overlap at $tstamp.") end
            end
        end
    end
end
