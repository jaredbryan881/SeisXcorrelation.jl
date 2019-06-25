using Geodesy
export rms, dist

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
