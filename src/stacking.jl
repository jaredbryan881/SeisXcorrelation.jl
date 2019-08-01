using SeisIO, JLD2, SeisNoise, Statistics, PlotlyJS
include("utils.jl")

"""
    threshold_stacking(data::CorrData, reference::CorrData, threshold::Float64)

Stack the windows in a CorrData object that exceed a correlation-coefficient threshold with respect to a reference.

# Arguments
- `data::CorrData,`    : Input CorrData matrix used to define the reference and each window
- `reference::CorrData`    : Reference cross-correlation function to use in computing the correlation coefficient
- `threshold::Float64,`    : Value of correlation-coefficient below which we zero windows

# Output
- `stackedData::CorrData,`    : Slectively stacked data
- `ccList::Array{Float64,1}`    : Correlation coefficient of each window with respect to the reference
"""
function selective_stacking(data::CorrData, reference::CorrData; threshold::Float64=0.0)
    # compute correlation coefficient for each window in data with respect to reference
    ccList = get_cc(data, reference)

    # find cross-correlations that fall below the cc threshold
    good_fit = findall(x->(x>=threshold), ccList)

    # copy data to extract well-fitting windows
    tempData = copy(data)
    tempData.corr = tempData.corr[:, good_fit]

    # linearly stack all data that exceeds the correlation-coefficient threshold
    stackedData = stack(tempData, allstack=true)

    return stackedData, ccList
end

"""
    get_cc(data::CorrData, reference::CorrData)

Compute the correlation coefficient for each window in a CorrData object with respect to a reference

# Arguments
- `data::CorrData`    : Input CorrData matrix used to define the reference and each window
- `reference::CorrData`    : Reference cross-correlation function to use in computing the correlation coefficient

# Output
- `ccList::Array{Float64,1}`    :  Correlation coefficient for each window with respect to the reference
"""
function get_cc(data::CorrData, reference::CorrData)
    # array of correlation coefficients for each constituent cross-correlation window
    ccList = Array{Float32,1}(undef,size(data.corr)[2])
    # iterate over each windowed cross-correlation and compute correlation-coefficient
    for j=1:size(data.corr)[2]
        ccList[j] = cor(reference.corr[:,1], data.corr[:,j])
    end

    return ccList
end
