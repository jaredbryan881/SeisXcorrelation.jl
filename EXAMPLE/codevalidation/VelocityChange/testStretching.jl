using SeisIO, Noise, JLD2, PlotlyJS

# This script tests the functionality of Stretching.jl from Noise.jl
# This assumes data generated via makeStretchData.jl

finame = "StretchedData"
valData = jldopen(finame*".jld2")

ref = valData["u0"]
cur = valData["u1"]
st = valData["st"]
dt = valData["dt"]
true_dvv = -0.01

# MWCS input parameters
InputDict = Dict( "freqmin"       => 0.001,
                  "freqmax"       => 9.9,
                  "fs"            => 20.0,
                  "mintimelag"    => 0.0,
                  "dtt_width"     => 30.0,
                  "dvmin"         => -0.015,
                  "dvmax"         => 0.015,
                  "ntrial"        => 100,
                  "window_length" => 30.0,
                  "window_step"   => 0)

time_axis = collect(0:1/InputDict["fs"]:300)

# allow windowed stretching analysis
if InputDict["window_step"] == 0
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
    dv, cc, cdp, eps, err, C, dtfiner = stretching(ref,
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
println("True: $true_dvv%")
println("Estimated: $(sum(dv_list)/N)%")
error = (sum(dv_list)/N) - true_dvv
println("Error: $error")
close(valData)
