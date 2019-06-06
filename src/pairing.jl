export sort_pairs
"""

    sort_pairs(pairs::AbstractArray)

Sort station pairs into their correlation category (auto-correlation, cross-correlation, cross-channel-correlation)

# Arguments
- `pairs::AbstractArray,`    : Array of station pairs, e.g. output of generate_pairs(file["info/stationlist"])

# Output
- `sorted_pairs::Dict{String, Array{String, N}}`    : dictionary containing station pairs grouped by correlation type

"""
function sort_pairs(pairs::AbstractArray)
    # dict containing each correlation category
    # preallocation of dict size assumes all stations have the same number of channels, which may not be true.
    sorted_pairs = Dict{String, Array{String}}("acorr"     => ["", ""],
                                               "xcorr"     => ["", ""],
                                               "xchancorr" => ["", ""])

    # fill dictionary based on detected correlation type
    for stnPair in pairs
        # same station, same channel
        if stnPair[1] == stnPair[2]
            sorted_pairs["acorr"] = hcat(sorted_pairs["acorr"], stnPair)
        # same station, different channel
        elseif (stnPair[1][end-3:end] != stnPair[2][end-3:end]) && (stnPair[1][1:end-3] == stnPair[2][1:end-3])
            sorted_pairs["xchancorr"] = hcat(sorted_pairs["xchancorr"], stnPair)
        # different station
        else
            sorted_pairs["xcorr"] = hcat(sorted_pairs["xcorr"], stnPair)
        end
    end

    # Remove ["", ""] used to initialized the dict
    sorted_pairs["acorr"]     = sorted_pairs["acorr"][:, 2:end]
    sorted_pairs["xcorr"]     = sorted_pairs["xcorr"][:, 2:end]
    sorted_pairs["xchancorr"] = sorted_pairs["xchancorr"][:, 2:end]

    return sorted_pairs
end
