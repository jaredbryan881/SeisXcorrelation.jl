"""
    save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)

    save CorrData structure to JLD2
"""
function save_CorrData2JLD2(foname::String, varname::String, CD::CorrData)
    file = jldopen(foname, "r+")
    file[varname] = CD
    JLD2.close(file)
end
