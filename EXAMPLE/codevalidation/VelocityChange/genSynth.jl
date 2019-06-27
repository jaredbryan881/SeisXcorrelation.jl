include("makeStretchData.jl")

using PlotlyJS

dvV = -0.01
dt = 0.05

foname="verificationData.jld2"
f=jldopen(foname, "a+")
dvVlist = collect(-0.01:0.0001:0.01)
noiselist = collect(0.0:0.001:0.1)
f["info/dvVlist"] = dvVlist
f["info/noiselist"] = noiselist
for dvV in dvVlist
    for noiselvl in noiselist
        # Example of damped sinusoid generation, stretching, and noise addition
        dampedSinParams = Dict( "A"    => 1.0,
                                "ω"    => 0.25,
                                "ϕ"    => 0.0,
                                "λ"    => 0.01,
                                "dt"   => 0.05,
                                "η"    => 150.0,
                                "t0"   => 0.0,
                                "npts" => 6001)

        signal1, t = generateSignal("dampedSinusoid", params=dampedSinParams)
        signal2, t = stretchData(signal1, dt, dvV, n=noiselvl*10)

        signal1 = addNoise(signal1, noiselvl)
        signal2 = addNoise(signal2, noiselvl)

        f["dampedSinusoid/$dvV.$noiselvl"] = [signal1, signal2]

        #=
        trace1 = scatter(;x=t, y=signal1, mode="lines")
        trace2 = scatter(;x=t, y=signal2, mode="lines")
        plots=[trace1, trace2]
        p=PlotlyJS.plot(plots)
        display(p)
        readline()
        =#

        # Example of ricker wavelet generation and convolution with random reflectivity
        # series, stretching, and noise addition.
        rickerParams = Dict( "f"     => 0.25,
                             "dt"    => 0.05,
                             "npr"   => 6001,
                             "npts"  => 6001)

        signal1, t = generateSignal("ricker", sparse=100, params=rickerParams)
        signal2, st = stretchData(signal1, dt, dvV, n=noiselvl*10)

        signal1 = addNoise(signal1, noiselvl)
        signal2 = addNoise(signal2, noiselvl)

        f["rickerConv/$dvV.$noiselvl"] = [signal1, signal2]

        #=
        trace1 = scatter(;x=t, y=signal1, mode="lines")
        trace2 = scatter(;x=t, y=signal2, mode="lines")
        plots=[trace1, trace2]
        p=PlotlyJS.plot(plots)
        display(p)
        readline()
        =#

    end
end
close(f)
