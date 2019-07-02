using JLD2
export read_JLD22Dict, save_CorrData2JLD2, save_Dict2JLD2, save_Array2JLD2

"""
    read_JLD22Dict(inFile, times::Array{Int64,1})

    read data structures from JLD2 into a list of dictionaries
"""
function read_JLD22Dict(inFile, tstamplist, times::Array{Int64,1})
    # create a list of dictionaries containing input data for each time step
    dsets = [Dict("$ts/$stn" => inFile["$ts/$stn"] for stn in keys(inFile[ts])) for ts in tstamplist[times]]
    return dsets
end

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
