using SeisIO, SeisNoise, JLD2, PlotlyJS

# This script tests the functionality of Stretching.jl from Noise.jl
# This assumes data generated via genSynth.jl

function testStretching(finame::String, foname::String, InputDict::Dict, windowlenlist::Array{Float64,1}, winsteplist::Array{Float64,1}, type::String)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    f_err = jldopen(foname, "a+")
    f_err["info/dvV"] = dvVlist
    f_err["info/noise"] = noiselist

    for dvV in dvVlist[1:5:end]
        println("dv/v: $dvV")
        for noiselvl in noiselist[:]
            for win_len in windowlenlist, win_step in windowsteplist
                signals = valData["$type/$dvV.$noiselvl"]
                (ref, cur) = signals

                time_axis = collect(InputDict["mintimelag"]:1/InputDict["fs"]:InputDict["mintimelag"]+InputDict["dtt_width"])

                # allow windowed stretching analysis
                if win_step == 0.0
                    minind = [1]
                    window_length_samples = length(time_axis)
                else
                    window_length_samples = convert(Int,win_len*InputDict["fs"])
                    window_step_samples = convert(Int, win_step*InputDict["fs"])
                    minind = 1:window_step_samples:length(ref) - window_length_samples
                end

                # number of windows
                N = length(minind)
                dv_list = zeros(N)
                err_list = zeros(N)

                # iterate over windows, report average dv/v
                for ii=1:N
                    window = collect(minind[ii]:minind[ii]+window_length_samples-1)
                    dv, cc, cdp, eps, err, C = stretching(ref,
                                                          cur,
                                                          time_axis,
                                                          window,
                                                          InputDict["freqmin"],
                                                          InputDict["freqmax"],
                                                          dvmin=InputDict["dvmin"],
                                                          dvmax=InputDict["dvmax"])
                    dv_list[ii] = dv
                    err_list[ii] = err
                end
                comp_dvV = sum(dv_list)/N
                comp_err = sum(err_list)/N

                f_err["$type/$dvV.$noiselvl"] = comp_err
            end
        end
    end
    close(f_err)
    close(valData)
end

# MWCS input parameters
InputDict_real = Dict( "freqmin"    => 0.1,
                       "freqmax"    => 9.9,
                       "fs"         => 20.0,
                       "mintimelag" => -100.0,
                       "dtt_width"  => 200.0,
                       "dvmin"      => -0.03,
                       "dvmax"      => 0.03,
                       "ntrial"     => 100 )

InputDict_synth = Dict( "freqmin"    => 0.1,
                        "freqmax"    => 9.9,
                        "fs"         => 20.0,
                        "mintimelag" => 0.0,
                        "dtt_width"  => 200.0,
                        "dvmin"      => -0.03,
                        "dvmax"      => 0.03,
                        "ntrial"     => 100 )

finame = "verificationData.jld2"
foname = "verificationDataError_stretching.jld2"
windowlenlist = [200.0]
windowsteplist = [0.0]
type = "rickerConv"

testStretching(finame, foname, InputDict_synth, windowlenlist, windowsteplist, type)
