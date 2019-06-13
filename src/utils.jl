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
    file = jldopen(foname, "r+")
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
