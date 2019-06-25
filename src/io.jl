using JLD2
export save_CorrData2JLD2, save_Dict2JLD2, save_Array2JLD2

"""
    save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)

    save CorrData structure to JLD2
"""
function save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)
    file = jldopen(foname, "a+")
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
