using PlotlyJS, Random, DSP, Dierckx, FileIO, JLD2

include("waves.jl")

# This script generates or reads a synthetic time series, stretches it by some
# value, and (optionally) adds noise to the signal.

"""

    generateSignal(type::String; params::Dict{String,Real}=Dict(), seed::Int64=66743)

Generate a signal of type "ricker" or "dampedSinusoid".

# Arguments
- `type::String,`    : Type of signal to generate: "ricker" or "dampedSinusoid"
- `params::Dict{String,Real}`    : parameters used to define the Ricker wavelet or damped sinusoid
- `seed::Int64`    : Seed for random number generator

# Output
- `u0`::Array{Float64,1}    : Computed signal
- `t`::Array{Float64,1}    : Time axis

"""
function generateSignal(type::String; params::Dict{String,Real}=Dict(), seed::Int64=66473)
    if type=="ricker"
        # unpack parameters
        f = params["f"]
        dt = params["dt"]
        npr = params["nptsRicker"]
        npts = params["npts"]

        # generate a ricker wavelet
        w, t = ricker(f=f, n=npr, dt=dt)

        # generate a random reflectivity series
        rng = MersenneTwister(seed)
        f = randn(rng, Float64, npts)

        # convolve reflectivity series with ricker wavelet
        u0 = conv(f,w)[1:npts]

    elseif type=="dampedSinusoid"
        # unpack parameters
        A = params["A"]
        ω = params["ω"]
        ϕ = params["ϕ"]
        λ = params["λ"]
        dt = params["dt"]
        npts = params["npts"]

        u0, t = dampedSinusoid(A=A, ω=ω, ϕ=ϕ, n=npts, dt=dt, λ=λ)
    end

    return u0, t
end


"""

    stretchData(u0::Array{Float64,1}, dt::Float64, dvV::Float64; n::Float64=0.0, seed::Int64=66743)

Linearly stretch a given signal, u0, by some factor given by a homogenous relative velocity change dvV.

# Arguments
- `u0::Array{Float64,1},`    : Input signal to be stretched
- `dt::Float64`    : Sampling interval
- `dvV::Float64`    : Stretching factor's causal relative velocity change
- `n::Float64`    : Noise level as a scalar of maximum time
- `seed::Int64`    : Seed for random number generator

# Output
- `u1`::Array{Float64,1}    : Stretched signal
- `st`::Array{Float64,1}    : Stretched time axis

"""
function stretchData(u0::Array{Float64,1}, dt::Float64, dvV::Float64; n::Float64=0.0, seed::Int64=66473)
    tvec = collect(( 0:length(u0)-1) .* dt)
    st = tvec .* dvV

    if (n != 0.0) st = addNoise(st, n) end

    tvec2 = tvec .+ st

    spl = Spline1D(tvec, u0)
    u1 = spl(tvec2)

    return u1, st
end

"""

    addNoise(signal::Array{Float64,1}, level::Float64; seed::Int64=66743)

Add noise to an array by some given mulitple of the maximum value.

# Arguments
- `signal::Array{Float64,1},`    : Input signal to add noise to
- `level::Float64`    : Multiple of maximum value of signal to set noise amplitude
- `seed::Int64`    : Seed for random number generator

# Output
- `signal`::Array{Float64,1}    : Signal with added noise

"""
function addNoise(signal::Array{Float64,1}, level::Float64; seed::Int64=66473)
    rng = MersenneTwister(seed)
    noise = randn(rng, Float64, length(signal))
    noise = level * maximum(abs.(signal)).*(noise./maximum(abs.(signal)))
    signal = signal .+ noise

    return signal
end

# Example of damped sinusoid generation, stretching, and noise addition
rickerParams = Dict( "f" => 2.5,
                     "dt" => 0.05,
                     "nptsRicker" => 41,
                     "npts" => 6001)

dampedSinParams = Dict( "A" => 1.0,
                        "ω" => 0.5,
                        "ϕ" => 0.0,
                        "λ" => 0.025,
                        "dt" => 0.05,
                        "npts" => 6001)

signal1, t = generateSignal("ricker", params=rickerParams)
signal2, st, dt = stretchData(signal1, 0.05, -0.1, n=0.05)

signal1 = addNoise(signal1, 0.001)
signal2 = addNoise(signal2, 0.001)

trace1 = scatter(;x=t, y=signal1, mode="lines")
trace2 = scatter(;x=t, y=signal2, mode="lines")
plots=[trace1, trace2]
p=PlotlyJS.plot(plots)
display(p)
readline()
