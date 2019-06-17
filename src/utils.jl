using Geodesy

"""
    save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)

    save CorrData structure to JLD2
"""
function save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)
    file = jldopen(foname, "r+")
    file[varname] = CD
    JLD2.close(file)
end

"""
    save_Dict2JLD2(foname::String, varname::String, D::Dict)

    save Dictionary to JLD2
"""
function save_Dict2JLD2(foname::String, varname::String, D::Dict)
    file = jldopen(foname, "r+")
    file[varname] = D
    JLD2.close(file)
end

"""
    save_Dict2JLD2(foname::String, varname::String, D::Dict)

    save Dictionary to JLD2
"""
function save_Array2JLD2(foname::String, varname::String, A::AbstractArray)
    file = jldopen(foname, "a+")
    file[varname] = A
    JLD2.close(file)
end

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
