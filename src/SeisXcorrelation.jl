__precompile__()
#module SeisXcorrelation

using SeisIO, SeisDownload, Noise, Printf, Dates, FFTW, JLD2, Distributed

include("pairing.jl")
include("fft.jl")
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

@everywhere function seisxcorrelation(tstamp::String, finame::String, foname::String, corrtype::Array{String,1}, corrorder::Int, maxtimelag::Real, freqmin::Real, freqmax::Real, fs::Real, cc_len::Int, cc_step::Int, IsAllComponents::Bool=false)
    # iterate over correlation categories
    for ct in corrtype
        # iterate over station pairs in the current correlation category
        for j = 1:length(sorted_pairs[ct][1, :])
            # get station names
            stn1 = sorted_pairs[ct][1, j]
            stn2 = sorted_pairs[ct][2, j]
            println("Corrrelating $stn1 with $stn2")

            # compute the xcorr on a pair of time series
            if corrorder == 1
                # read station SeisChannels into SeisData before FFT
                S1 = SeisData(data["$tstamp/$stn1"])
                if ct=="acorr" S2=S1 else S2=SeisData(data["$tstamp/$stn2"]) end # S2 is a ref to S1 if "acorr"

                # compute FFT using Noise.jl -- returns type FFTData
                FFT1 = compute_fft(S1, freqmin, freqmax, fs, cc_step, cc_len)
                if ct=="acorr" FFT2=FFT1 else FFT2=compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len) end # FFT2 is a ref to FFT1 if "acorr"

                # compute correlation using Noise.jl -- returns type CorrData
                xcorr = compute_cc(FFT1, FFT2, maxtimelag)

            # compute the xcorr on a pair of xcorrelations
            elseif corrorder == 2
                # read CorrData
                C1 = data["$tstamp/$stn1"]
                C2 = data["$tstamp/$stn2"]

                # compute FFT on CorrData direcly
                FFT1 = compute_fft_c3(C1, freqmin, freqmax, fs, cc_step, cc_len)

                println("Higher order correlation methods are not yet implemented.\nExiting code.")
                exit()

            # compute xcorr on a pair of xcorrelation codas
            elseif corrorder == 3
                println("Higher order correlation methods are not yet implemented.\nExiting code.")
                exit()
            end

            # save data after each cross correlation
            varname = "$tstamp/$stn1.$stn2"
            save_CorrData2JLD2(foname, varname, xcorr)
        end
    end
    return 0
end

"""
    save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)

    save CorrData structure to JLD2
"""
function save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)
    file = jldopen(foname, "r+")
    file[varname] = CD
    JLD2.close(file)
end

#end
