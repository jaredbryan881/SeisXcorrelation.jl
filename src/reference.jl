using SeisIO, JLD2, Noise
include("utils.jl")

function compute_reference_xcorr(finame::String, foname::String, xcorrSize::Int)
    f = jldopen(finame)

    corrstationlist = f["info/corrstationlist"] # names of station pairs
    tslist = f["info/timestamplist"] # time stamps

    jldopen(foname, "w") do file
        file["info/timestamplist"]   = tslist;
        file["info/corrstationlist"] = corrstationlist;
    end

    # iterate over correlation type
    for ct in keys(corrstationlist)
        # iterate over station pairs
        for stnum=1:size(corrstationlist[ct])[2]
            # get name of station pair
            stn1 = corrstationlist[ct][1, stnum]
            stn2 = corrstationlist[ct][2, stnum]
            varname = "$stn1" .* "." .* "$stn2"
            println(varname)

            # initialize array to hold stacked data
            stackData = zeros(xcorrSize)
            # num of summed xcorrs
            numxcorr = 0
            # compute full stack for this station pair
            for tstamp in tslist
                if stn1 in f["$tstamp/ccerrors"]["DimensionMismatch"] || stn2 in f["$tstamp/ccerrors"]["DimensionMismatch"]
                    println("Passing")
                    continue
                elseif stn1 in f["$tstamp/ccerrors"]["DataUnavailable"] || stn2 in f["$tstamp/ccerrors"]["DataUnavailable"]
                    println("Passing")
                    continue
                end
                # read data and stack
                data = f["$tstamp/$varname"].corr
                stackData += sum(data, dims=2)
                numxcorr += size(data)[2]
            end
            # compute average
            stackData /= numxcorr
            # save reference xcorr
            save_Array2JLD2(foname, "$varname", stackData)
        end
    end
end

compute_reference_xcorr("../EXAMPLE/xcorr_BP/outputData/BPnetworkxcorr_weq.jld2", "reference_xcorr.jld2", 4001)
