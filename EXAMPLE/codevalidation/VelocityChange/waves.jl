export ricker, dampedSinusoid

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
function ricker(;f::Float64=20.0, n::Int64=100, dt::Float64=0.001)
    # Create the wavelet and shift in time if needed
    T = dt*(n-1);
    t = collect(0:dt:T);
    t0 = 1/f;
    tau = t.-t0;

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
function dampedSinusoid(;A::Float64=1.0, ω::Float64=1.0, ϕ::Float64=0.0, η::Float64=0.0, n::Int64=100, dt::Float64=0.001, λ::Float64=2.0)
    T = dt * (n-1)
    t = collect(0:dt:T)

    dampedSin = A*exp.(-λ*abs.(t.-η)) .* cos.(ω*t .+ ϕ)

    return (dampedSin, t)
end
