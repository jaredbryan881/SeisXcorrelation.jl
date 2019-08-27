using SeisIO, SeisNoise, Geodesy, SeisIO
export rms, dist, normalize, normalize!

"""
    rms(A::AbstractArray, B::AbstractArray)

Compute the root mean squared error between two arrays
"""
function rms(A::AbstractArray, B::AbstractArray)
    return sqrt(sum((B.-A).^2)/length(A))
end

"""
    dist(loc1::GeoLoc, loc2::GeoLoc)

Compute the distance (m) between two points (lat, lon, elev)
This assumes the distance is not large, such that the curvature of the earth is negligible.
"""
function dist(loc1::GeoLoc, loc2::GeoLoc)
    loc1_lla = LLA(loc1.lat, loc1.lon, loc1.el)
    loc2_lla = LLA(loc2.lat, loc2.lon, loc2.el)

    return distance(loc1_lla, loc2_lla)
end

"""
    window_read(len::Int64, max_len::Int64; startInd::Int64=2)

Define subarrays of a given length for data reading.
"""
function window_read(len::Int64, max_len::Int64)
    fullArr = collect(1:len)
    # overflow
    rem = length(fullArr) % max_len
    # number of full sized subarrays
    nWin = div(length(fullArr), max_len)

    # define subarrays
    windows = [fullArr[1+(i*max_len):(i+1)*max_len] for i=0:nWin-1]
    if rem != 0 push!(windows, collect((max_len*nWin)+1:len)) end

    return windows
end

"""
    normalize!(x::AbstractArray; mag::Float64=1.0)

Normalize an array and scale to a given magnitude.
"""
function normalize!(x::AbstractArray; mag::Float64=1.0)
    x ./= maximum(abs.(x)) .* mag
    return nothing
end
normalize(x::AbstractArray; mag::Float64=1.0) = (U = deepcopy(x); normalize!(U, mag=mag); return U)

"""
    Base.copy(x::T)

Universal copy function used to copy composite types such as structs.
"""
Base.copy(x::T) where T = T([deepcopy(getfield(x, k)) for k ∈ fieldnames(T)]...)
