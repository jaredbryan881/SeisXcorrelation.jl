export generateSignal, stretchData, addNoise, addNoise!

using Random, DSP, Dierckx, FileIO, JLD2, SeisIO, FFTW

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
function generateSignal(type::String, params::Dict{String,Real}; sparse::Int64=0, seed::Int64=66473, stretchSource::Float64=0.0)
    if type=="ricker"
        # unpack parameters
        f    = params["f"]    # peak frequency
        dt   = params["dt"]   # sampling interval
        npr  = params["npr"]  # number of points in the Ricker wavelet
        npts = params["npts"] # number of points in the generated signal
        m    = params["m"]    # power for reflectivity series

        # generate a ricker wavelet
        (w, t) = ricker(f=f, n=npr, dt=dt)

        if stretchSource != 0.0
            spectrum = rfft(w)
            stAmpSpec, st = stretchData(real(spectrum), dt, stretchSource)
            spectrum = stAmpSpec .+ imag(spectrum).*im
            w = irfft(spectrum, length(w))
        end

        # generate a random reflectivity series
        rng = MersenneTwister(seed)
        f = randn(rng, Float64, npts)

        if sparse!=0
            r=collect(1:npts)
            f[r .% sparse .!= 0] .= 0
        end
        # raise output from randn to an integral power. See <<Numerical Methods of
        # Explorational Seismology>> page 77.
        f = sign.(f) .* abs.(f).^m

        # convolve reflectivity series with ricker wavelet
        u0 = conv(f,w)[convert(Int64, npts-floor(npts/2)):convert(Int64, npts+floor(npts/2))]

    elseif type=="dampedSinusoid"
        # unpack parameters
        A    = params["A"]    # amplitude
        ω    = params["ω"]    # angular frequency of the sinusoid
        ϕ    = params["ϕ"]    # phase shift for the sinusoid
        η    = params["η"]    # shift in function maximum
        λ    = params["λ"]    # decay constant for the exponential decay
        dt   = params["dt"]   # sampling interval
        t0   = params["t0"]   # start time
        npts = params["npts"] # number of points in the generated signal

        u0, t = dampedSinusoid(A=A, ω=ω, ϕ=ϕ, η=η, n=npts, dt=dt, t0=t0, λ=λ)

    elseif type=="sinc"
        # unpack parameters
        A    = params["A"]
        ω    = params["ω"]    # angular frequency of the sinusoid
        ϕ    = params["ϕ"]    # phase shift for the sinusoid
        dt   = params["dt"]   # sampling interval
        t0   = params["t0"]   # start time
        npts = params["npts"] # number of points in the generated signal

        u0, t = sinc(A=A, ω=ω, ϕ=ϕ, dt=dt, t0=t0, n=npts)

    elseif type=="chirp"
        # unpack parameters
        c     = params["c"]     # phase velocity
        tp    = params["tp"]    # period
        mintp = params["mintp"] # minimum period
        maxtp = params["maxtp"] # maximum period
        dist  = params["dist"]  # distance
        npts  = params["npts"]  # number of points in the generated signal
        dt    = params["dt"]    # sampling interval
        t0    = params["t0"]    # start time

        u0, t = chirp(c=c, tp=tp, mintp=mintp, maxtp=maxtp, dist=dist, n=npts, dt=dt, t0=t0)

    elseif type=="spectSynth"
        # unpack parameters
        A    = params["A"]
        dt   = params["dt"]
        t0   = params["t0"]
        npts = params["npts"]

        u0, t = spectSynth(A=A, dt=dt, t0=t0, n=npts)
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
- `tvec`::Array{Float64,1}    : Unstretched time axis

"""
function stretchData(u0::Union{Array{Float32,1},Array{Float64,1}}, dt::Float64, dvV::Float64; starttime::Float64=0.0, stloc::Float64=0.0, n::Float64=0.0, seed::Int64=66473)
    tvec = collect((0:length(u0)-1) .* dt) .+ starttime
    st = (tvec.-stloc) .* dvV

    if (n != 0.0) addNoise!(st, n, seed=seed) end

    tvec2 = tvec .+ st

    # interpolation
    spl = Spline1D(tvec, u0)
    u1 = spl(tvec2)

    return u1, tvec
end

"""

    addNoise!(signal::Array{Float64,1}, level::Float64; seed::Int64=66743)

Add noise to an array by some given mulitple of the maximum value.

# Arguments
- `signal::Array{Float64,1},`    : Input signal to add noise to
- `level::Float64`    : Multiple of maximum value of signal to set noise amplitude
- `seed::Int64`    : Seed for random number generator

# Output
- `signal`::Array{Float64,1}    : Signal with added noise

"""
function addNoise!(signal::Union{Array{Float32,1},Array{Float64,1}}, level::Float64; percent::Bool=true, seed::Int64=66473, freqmin::Float64=0.1, freqmax::Float64=9.99, fs=20.0, corners::Int64=4, zerophase::Bool=false)
    rng = MersenneTwister(seed)
    noise = randn(rng, Float64, length(signal))
    bandpass!(noise, freqmin, freqmax, fs, corners=corners, zerophase=zerophase)
    noise ./= maximum(abs.(noise)) # normalize noise
    noise .*= level # scale noise to given level
    if percent noise .*= maximum(abs.(signal)) end
    signal .+= noise

    return nothing
end
addNoise(signal::Union{Array{Float32,1},Array{Float64,1}}, level::Float64;
         seed::Int64=66473, freqmin::Float64=0.1, freqmax::Float64=9.99, fs=20.0,
         corners::Int64=4, zerophase::Bool=false) = (U = deepcopy(signal);
         addNoise!(U, level, seed=seed, freqmin=freqmin, freqmax=freqmax, fs=fs,
         corners=corners, zerophase=zerophase); return U)
