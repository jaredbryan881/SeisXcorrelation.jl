export generateSignal, stretchData, addNoise

using Random, DSP, Dierckx, FileIO, JLD2

include("waves.jl")

# This script generates or reads a synthetic time series, stretches it by some
# value, and (optionally) adds noise to the signal.

"""

    generateSignal(type::String; params::Dict{String,Real}=Dict(), sparse::Int64=0, seed::Int64=66743)

Generate a signal of type "ricker" or "dampedSinusoid".

# Arguments
- `type::String,`    : Type of signal to generate: "ricker" or "dampedSinusoid"
- `params::Dict{String,Real}`    : parameters used to define the Ricker wavelet or damped sinusoid
- `sparse::Int64`    : Interval over which to keep elements in array, e.g. sparse=4 corresponds to keeping every 4th element
- `seed::Int64`    : Seed for random number generator

# Output
- `u0`::Array{Float64,1}    : Computed signal
- `t`::Array{Float64,1}    : Time axis

"""
function generateSignal(type::String; params::Dict{String,Real}=Dict(), sparse::Int64=0, seed::Int64=66473)
    if type=="ricker"
        # unpack parameters
        f    = params["f"]    # peak frequency
        dt   = params["dt"]   # sampling interval
        npr  = params["npr"]  # number of points in the Ricker wavelet
        npts = params["npts"] # number of points in the generated signal

        # generate a ricker wavelet
        (w, t) = ricker(f=f, n=npr, dt=dt)

        # generate a random reflectivity series
        rng = MersenneTwister(seed)
        f = randn(rng, Float64, npts)
        if sparse!=0
            r=collect(1:npts)
            f[r .% sparse .!= 0] .= 0
        end

        # convolve reflectivity series with ricker wavelet
        u0 = conv(f,w)[1:npts]

    elseif type=="dampedSinusoid"
        # unpack parameters
        A    = params["A"]    # amplitude
        ω    = params["ω"]    # angular frequency of the sinusoid
        ϕ    = params["ϕ"]    # phase shift for the sinusoid
        η    = params["η"]    # shift in function maximum
        λ    = params["λ"]    # decay constant for the exponential decay
        dt   = params["dt"]   # sampling interval
        t0   = params["t0"]   # leftmost time
        npts = params["npts"] # number of points in the generated signal

        u0, t = dampedSinusoid(A=A, ω=ω, ϕ=ϕ, η=η, n=npts, dt=dt, t0=t0, λ=λ)
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
function stretchData(u0::Array{Float64,1}, dt::Float64, dvV::Float64; starttime::Float64=0.0, stloc::Float64=0.0, n::Float64=0.0, seed::Int64=66473)
    tvec = collect(( 0:length(u0)-1) .* dt) .+ starttime
    st = (tvec.-stloc) .* dvV

    if (n != 0.0) st = addNoise(st, n) end

    tvec2 = tvec .+ st

    # interpolation
    spl = Spline1D(tvec, u0)
    u1 = spl(tvec2)

    return u1, tvec
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
