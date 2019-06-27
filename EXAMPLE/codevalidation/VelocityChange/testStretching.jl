using SeisIO, Noise, JLD2, PlotlyJS

# This script tests the functionality of Stretching.jl from Noise.jl
# This assumes data generated via genSynth.jl

function testStretching(finame::String, windowlenlist::Array{Float64,1}, winsteplist::Array{Float64,1}, type::String)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    for dvV in dvVlist[1:10:end]
        for noiselvl in noiselist[1:1]
            for win_len in windowlenlist, win_step in windowsteplist
                signals = valData["$type/$dvV.$noiselvl"]
                ref = signals[1]
                cur = signals[2]

                # MWCS input parameters
                InputDict = Dict( "freqmin"       => 0.1,
                                  "freqmax"       => 9.9,
                                  "fs"            => 20.0,
                                  "mintimelag"    => 0.0,
                                  "dtt_width"     => 300.0,
                                  "dvmin"         => -0.015,
                                  "dvmax"         => 0.015,
                                  "ntrial"        => 100,
                                  "window_length" => win_len,
                                  "window_step"   => win_step)

                time_axis = collect(0:1/InputDict["fs"]:300)

                # allow windowed stretching analysis
                if InputDict["window_step"] == 0.0
                    minind = [1]
                    window_length_samples = length(time_axis)
                else
                    window_length_samples = convert(Int,InputDict["window_length"] * InputDict["fs"])
                    window_step_samples = convert(Int,InputDict["window_step"] * InputDict["fs"])
                    minind = 1:window_step_samples:length(ref) - window_length_samples
                end

                N = length(minind)
                dv_list = zeros(N)

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
                end

                # Report performance
                println("True: $(dvV*100)%")
                println("Estimated: $(sum(dv_list)/N)%")
                error = (sum(dv_list)/N) - (dvV*100)
                println("Error: $error%")
                println("")
            end
        end
    end
    close(valData)
end

finame = "verificationData.jld2"
windowlenlist = [300.0]
windowsteplist = [0.0]
type = "rickerConv"
testStretching(finame, windowlenlist, windowsteplist, type)
