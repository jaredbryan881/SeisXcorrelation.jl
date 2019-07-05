using SeisIO, SeisNoise, DTWDT, PlotlyJS, JLD2
import LinearAlgebra.pinv

# This is a script to test the functionality of the DTW functions of DTWDT.jl
# This assumes synthetic data generated via genSynth.jl

function testDTW(finame::String, InputDict::Dict, type::String)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    for dvV in dvVlist[1:10:end]
        for noiselvl in noiselist[1:10:end]
            signals = valData["$type/$dvV.$noiselvl"]
            (u1, u2) = signals

            maxLag = InputDict["maxLag"]
            b = InputDict["b"]
            st = InputDict["st"]
            dt = InputDict["dt"]
            direction = InputDict["direction"]
            mintime = InputDict["mintime"]

            npts = length(u1)
            tvec = (collect(0:npts-1) .* dt) .+ mintime

            (stbarTime, stbar, dist, error) = dtwdt(u1, u2, dt, maxLag=maxLag, b=b, direction=direction)

            A = ones(npts,2)
            A[:,1] = Array(1:npts) ./ npts
            coef = pinv(A) * stbarTime
            (a, m) = coef

            # Report results
            println("Noise level: $(noiselvl*100)%")
            println("True dvV: $(dvV*100)%")
            println("Regressed dvV: $(-m)%")
            println("Distance: $error")

            plot_dtt(tvec, stbarTime, -m/100, dvV)
        end
    end
end

function plot_dtt(tvec::Array{Float64,1}, stbarTime::Array{Float64,1}, m::Float64, dvV::Float64)
    trace1 = PlotlyJS.scatter(;x=tvec, y=stbarTime, name="Predicted dt")
    trace2 = PlotlyJS.scatter(;x=tvec, y=tvec .* m, name="Regression")
    trace3 = PlotlyJS.scatter(;x=tvec, y=tvec .* dvV, name="True dt")
    p=PlotlyJS.plot([trace1, trace2, trace3])
    display(p)
    readline()
end

# MWCS input parameters
InputDict = Dict( "maxLag"    => 150,
                  "b"         => 1,
                  "st"        => 0,
                  "dt"        => 0.05,
                  "direction" => 0,
                  "mintime"   => -100.0)

finame = "verificationData.jld2"
type = "realData"

testDTW(finame, InputDict, type)
