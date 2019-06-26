include("makeStretchData.jl")

using PlotlyJS

dvV = -0.1
dt = 0.05
# Example of damped sinusoid generation, stretching, and noise addition
dampedSinParams = Dict( "A"    => 1.0,
                        "ω"    => 0.25,
                        "ϕ"    => 0.0,
                        "λ"    => 0.025,
                        "dt"   => 0.05,
                        "η"    => 50.0,
                        "npts" => 6001)

signal1, t = generateSignal("dampedSinusoid", params=dampedSinParams)
signal2, t = stretchData(signal1, dt, dvV, n=0.1)

signal1 = addNoise(signal1, 0.001)
signal2 = addNoise(signal2, 0.001)

trace1 = scatter(;x=t, y=signal1, mode="lines")
trace2 = scatter(;x=t, y=signal2, mode="lines")
plots=[trace1, trace2]
p=PlotlyJS.plot(plots)
display(p)
readline()

# Example of ricker wavelet generation and convolution with random reflectivity
# series, stretching, and noise addition.
rickerParams = Dict( "f"     => 2.5,
                     "dt"    => 0.05,
                     "npr"   => 41,
                     "npts"  => 6001)

signal1, t = generateSignal("ricker", params=rickerParams)
signal2, st = stretchData(signal1, dt, dvV, n=0.05)

signal1 = addNoise(signal1, 0.001)
signal2 = addNoise(signal2, 0.001)

trace1 = scatter(;x=t, y=signal1, mode="lines")
trace2 = scatter(;x=t, y=signal2, mode="lines")
plots=[trace1, trace2]
p=PlotlyJS.plot(plots)
display(p)
readline()
