using JLD2
export read_JLD22Dict, save_CorrData2JLD2, save_Dict2JLD2, save_Array2JLD2

"""
    read_JLD22Dict(inFile::Array{JLD2.JLDFile,1}, tstamplist::AbstractArray, times::Array{Int64,1})

    read data structures from JLD2 into a list of dictionaries
"""
function read_JLD22Dict(inFile::Array{JLD2.JLDFile{JLD2.MmapIO},1}, tstamplist::AbstractArray)
    # create a list of dictionaries containing input data for each time step
    # this assumes inFile is in the same order as tstamplist
    dsets = [Dict("$ts/$stn" => inFile[i]["$ts/$stn"] for stn in keys(inFile[i][ts])) for (i, ts) in enumerate(tstamplist)]
    return dsets
end

"""
    read_JLD22Dict(inFile, tstamplist::AbstractArray, times::Array{Int64,1})

    read data structures from JLD2 into a list of dictionaries
"""
function read_JLD22Dict(inFile::JLD2.JLDFile, tstamplist::AbstractArray)
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
