using PlotlyJS, Random, DSP, Dierckx, FileIO, JLD2

include("waves.jl")

# This script generates or reads a synthetic time series, and stretches it by some value

function generateSignal(type::String, f::Float64, npts::Int64, dt::Float64; nptsRicker::Int64=41, seed::Int64=66473)
    if type=="ricker"
        # generate a ricker wavelet
        w, t = ricker(f=f, n=nptsRicker, dt=dt)

        # generate a random reflectivity series
        rng = MersenneTwister(seed)
        f = randn(rng, Float64, npts)

        # convolve reflectivity series with ricker wavelet
        u0 = conv(f,w)[1:npts]

    elseif type=="dampedSinusoid"
        A=1.0
        ϕ=0.0
        λ=0.025
        u0, t = dampedSinusoid(A=A, ω=f, ϕ=ϕ, n=npts, dt=dt, λ=λ)
    end

    return u0, t
end

function stretchData(u0::Array{Float64,1}, dt::Float64, dvV::Float64; n::Float64=0.0, seed::Int64=66473)
    tvec = collect(( 0:length(u0)-1) .* dt)
    st = tvec .* dvV

    if (n != 0.0) st = addNoise(st, n) end

    tvec2 = tvec .+ st

    spl = Spline1D(tvec, u0)
    u1 = spl(tvec2)

    return u1, st, dt
end

function addNoise(signal::Array{Float64,1}, level::Float64 ;seed::Int64=66473)
    rng = MersenneTwister(seed)
    noise = randn(rng, Float64, length(signal))
    noise = level * maximum(abs.(signal)).*(noise./maximum(abs.(signal)))
    signal = signal .+ noise

    return signal
end

# Example of damped sinusoid generation, stretching, and noise addition
dt = 0.05
f = 0.5
npts = 6001
dvV = - 0.05

signal1, t = generateSignal("dampedSinusoid", f, npts, dt)
signal2, st, dt = stretchData(signal1, dt, dvV, n=0.05)

signal1 = addNoise(signal1, 0.01)
signal2 = addNoise(signal2, 0.01)

trace1 = scatter(;x=t, y=signal1, mode="lines")
trace2 = scatter(;x=t, y=signal2, mode="lines")
plots=[trace1, trace2]
p=PlotlyJS.plot(plots)
display(p)
readline()

@save foname dt u0 u1 st
