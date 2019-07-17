export ricker, dampedSinusoid, sinc

"""

    ricker(;f::Float64=20.0, n::Int64=100, dt::Float64=0.001, nargout::Bool=false)

Generate a Ricker wavelet.

# Arguments
- `f::Float64,`    : Peak frequnecy
- `n::Int64`    : Number of points
- `dt::Float64`    : Sampling interval

# Output
- `s`::Array{Float64,1}    : Computed Ricker wavelet
- `t`::Array{Float64,1}    : Time axis

"""
function ricker(;f::Float64=20.0, n::Int64=41, dt::Float64=0.05)
    # Create the wavelet
    t = timeAxis(dt, 0., n)
    tau = t .- 1/f
    s = (1.0 .- tau.*tau.*f.^2*π.^2).*exp.(-tau.^2*π.^2*f.^2);

    return (s, t)
end

"""

    dampedSinusoid(;A::Float64=1.0, ω::Float64=1.0, ϕ::Float64=0.0, n::Int64=100, dt::Float64=0.001, λ::Float64=2.0)

Generate a damped sinusoid.

# Arguments
- `A::Float64`    : Amplitude of the damped sinusoid
- `n::Int64`    : Number of points
- `dt::Float64`    : Sampling interval
- `ω::Float64`    : Frequency
- `λ::Float64`    : Decay constant
- `ϕ::Float64`    : Phase angle at t=0
- `η::Float64`    : Shift in function maximum

# Output
- `dampedSin::Array{Float64,1}`    : Computed damped sinusoid
- `t::Array{Float64,1}`    : Time axis

"""
function dampedSinusoid(;A::Float64=1.0, ω::Float64=1.0, ϕ::Float64=0.0, η::Float64=0.0, n::Int64=100, dt::Float64=0.001, t0::Float64=0.0, λ::Float64=2.0)
    t = timeAxis(dt, t0, n)
    dampedSin = A*exp.(-λ*abs.(t.-η)) .* cos.(ω*t .+ ϕ)

    return (dampedSin, t)
end

"""

    sinc(;A::Float64=1.0, ω::Float64=1.0, ϕ::Float64=0.0, n::Int64=100, dt::Float64=0.001, t0=0.0)

Generate a sinc (cardinal-sin) function

# Arguments
- `A::Float64`    : Amplitude of the damped sinusoid
- `n::Int64`    : Number of points
- `dt::Float64`    : Sampling interval
- `ω::Float64`    : Frequency
- `ϕ::Float64`    : Phase angle at t=0

# Output
- `s::Array{Float64,1}`    : Computed sinc
- `t::Array{Float64,1}`    : Time axis

"""
function sinc(;A::Float64=1.0, ω::Float64=1.0, ϕ::Float64=0.0, n::Int64=100, dt::Float64=0.001, t0::Float64=0.0)
    t = timeAxis(dt, t0, n)
    s = A.*Base.sinc.(ω.*t.+ϕ)

    return (s, t)
end

"""

    timeAxis(dt::Float64, t0::Float64, n::Float64)

Create a vector of times given sample frequency, number of points, and start time.

# Arguments
- `dt::Float64`    : Sample interval
- `t0::Int64`    : Start time
- `n::Int64`    : Number of points in the time vector

# Output
- `t::Array{Float64,1}`    : Time axis

"""
function timeAxis(dt::Float64, t0::Float64, n::Int64)
    T = dt * (n-1)
    t = collect(0:dt:T) .+ t0
    return t
end
