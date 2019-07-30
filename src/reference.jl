using SeisIO, JLD2, SeisNoise, Statistics, PlotlyJS
include("utils.jl")
include("partition.jl")
include("io.jl")

"""
compute_reference_xcorr(finame::String, foname::String)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `basefiname::String,`    : Input base file name e.g. "./inputData/BPnetwork"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair

"""
function compute_reference_xcorr(basefiname::String, foname::String; phase_smoothing::Float64=0.)
    # input file holds metadata (stationlist, timestamplist, etc.)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"]
    close(f) # base xcorr file

    # hold reference xcorrs in memory and write all at once
    ref_dict = Dict()

    # iterate over timestamps
    for (t, tstamp) in enumerate(tslist)
        # read unstacked xcorrs for each time stamp
        f_cur = jldopen(basefiname*".$tstamp.jld2")
        grp = f_cur[tstamp] # xcorrs
        println("stacking xcorrs at $tstamp")

        # iterate over station pairs
        for pair in keys(grp)
            # load xcorr and reverse if necessary
            xcorr = try grp[pair] catch; continue end

            # stack xcorrs over DL_time_unit (over time length of CorrData object)
            stack!(xcorr, allstack=true, phase_smoothing=phase_smoothing)

            # stack xcorrs if they have a key, assign key if not
            if haskey(ref_dict, pair)
                ref_dict[pair].corr .+= xcorr.corr
            else
                ref_dict[pair] = xcorr
            end
        end
        close(f_cur) # current xcorr file
    end

    # write output file of reference xcorrs
    f_ref = jldopen(foname, "a+")
    for pair in keys(ref_dict)
        f_ref[pair] = ref_dict[pair]
    end
    close(f_ref) # reference xcorr file
end

"""
convergence(finame::String, foname::String; metric::String="snr", skipoverlap::Bool=false)

Compute the RMS error or SNR for progressively longer stacks of cross-correlations and (for RMS) a reference cross-correlation

# Arguments
- `finame::String,`    : Input file name for unstacked cross-correlations e.g. "./inputData/BPnetwork.jld2"
- `refname::String,`    : Input file name for reference cross-corrleations e.g. "./inputData/BPnetwork.jld2"
- `foname::String,`    : Output file name for fully stacked xcorrs e.g. "referenceXcorr.jld2"
- `metric::String`    : Method of tracking convergence: "snr" for peak signal-to-noise ratio, "rms" for root mean square misfit, or "cc" for correlation coefficient
- `skipoverlap::Bool`    : Whether or not to skip overlap in cross-correlation windows. If cross-correlations were computed
                           with no overlap, set this to false to ensure no time is skipped.
- `slice::Union{Bool, Float64, Array{Float64,1}}`    : whether to slice the cross-correlations before computing convergence. If so, a single Float64 specifies the
                                                       lag at which to slice to keep just the ballistic wave and an array of Float64 specifies the start and end lag
                                                       to keep only the coda.

# Output
- `foname.jld2`    : contains arrays of RMS errors for progressively longer stacks of cross-correlations

"""
function convergence(basefiname::String, refname::String, foname::String; ntimes::Int64=-1, metric::String="cc", skipoverlap::Bool=false, slice::Union{Bool, Float64, Array{Float64,1}}=false)
    # input file contains metadata (stationlist, timestamplist)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"] # time stamps
    close(f) # input file

    # read reference xcorrs
    ref_f = jldopen(refname)

    # dictionary to keep track of convergence vectors for each station pair
    conv = Dict()
    # dictionary to keep track of stacked data
    stackData = Dict{String, Array{Float64,2}}()

    if ntimes==-1
        ntimes=length(tslist)
    end

    # iterate over timestamps
    for tstamp in tslist[1:ntimes]
        println("Processing $tstamp")
        # read xcorrs for each time step
        f_cur = jldopen(basefiname*".$tstamp.jld2")

        # iterate over station pairs
        for stnpair in keys(f_cur[tstamp])
            # load unstacked cross-correlations
            data = try f_cur["$tstamp/$stnpair"] catch; continue end

            # load reference cross-correlation for the current station pair
            reference = ref_f[stnpair]

            # optionally slice reference to keep only coda or only ballistic
            if typeof(slice) != Bool
                if typeof(slice)==Float64
                    # time vector
                    tvec = -reference.maxlag:1/reference.fs:reference.maxlag
                    # find all times within [-slice,slice] time
                    t_inds = findall(x->(x>=-slice && x<=slice), tvec)
                    # slice reference (maintain 2d array)
                    ref = reference.corr[t_inds, :]
                elseif typeof(slice)==Array{Float64,1}
                    # convert startlag/endlag[s] to startlag/windowlength[samples]
                    win_len = Int(diff(slice)[1] * reference.fs)
                    startlag = Int(slice[1] * reference.fs)
                    # partition reference to keep only the coda
                    ref = partition(reference.corr, startlag, win_len)
                else
                    println("Please choose an allowable slicing operation. Exiting.")
                    exit()
                end
            else
                # default to full reference cross-correlation
                ref = reference.corr
            end

            # iterate over windowed cross-correlations
            for i=1:size(data.corr)[2]
                if skipoverlap i=2*i-1 end # assumes cc_step is 1/2 cc_len

                # slice current cross-correlation if necessary
                if typeof(slice) != Bool
                    if typeof(slice)==Float64
                        # keep only ballistic -- 1D array
                        data_cur = data.corr[t_inds, i]
                    elseif typeof(slice)==Array{Float64,1}
                        # keep only coda -- 2D array (this will change implementation of convergence calculation)
                        data_cur = partition(data.corr[:, i:i], startlag, win_len)
                    end
                else
                    # default to full windowed cross-correlation -- 1D array
                    data_cur = data.corr[:, i]
                end

                # build input to convergence functions
                # here we either create a progressively longer stack of cross-correlation
                # functions, or concatenate cross-correlation functions for snr
                if haskey(stackData, stnpair) && (metric=="rms" || metric=="cc")
                    # keep running sum of cross-correlations, stack if key exists, assign if not
                    stackData[stnpair] .+= data_cur
                elseif haskey(stackData, stnpair) && metric=="snr"
                    # create matrix of xcorrs
                    stackData[stnpair] = hcat(stackData[stnpair], data_cur)
                else
                    # reshape into 2D array (necessary for snr input)
                    if typeof(slice) == Array{Float64,1}
                        stackData[stnpair] = data_cur
                    else
                        stackData[stnpair] = reshape(data_cur, length(data_cur), 1)
                    end
                end

                # compute convergence data point and append to vector
                if metric=="rms"
                    # compute RMS error of the running sum with the reference
                    if haskey(conv, stnpair)
                        push!(conv[stnpair], rms(normalize(ref[:,1]), normalize(stackData[stnpair][:,1])))
                    else
                        conv[stnpair] = [rms(normalize(ref[:,1]), normalize(stackData[stnpair][:,1]))]
                    end
                elseif metric=="cc"
                    # compute Pearson correlation-coefficient between current and reference cross-correlation
                    if haskey(conv, stnpair)
                        if typeof(slice)==Array{Float64,1}
                            # compute correlation-coefficient for both positive and negative coda
                            # store in 2d array
                            neg_coda_cor = cor(ref[:, 1], stackData[stnpair][:, 1])
                            pos_coda_cor = cor(ref[:, 2], stackData[stnpair][:, 2])
                            push!(conv[stnpair], [neg_coda_cor, pos_coda_cor])
                        else
                            push!(conv[stnpair], cor(ref[:,1], stackData[stnpair][:,1]))
                        end
                    else
                        if typeof(slice)==Array{Float64,1}
                            # compute correlation-coefficient for both positive and negative coda
                            # store in 2d array
                            neg_coda_cor = cor(ref[:, 1], stackData[stnpair][:, 1])
                            pos_coda_cor = cor(ref[:, 2], stackData[stnpair][:, 2])
                            conv[stnpair] = [[neg_coda_cor, pos_coda_cor]]
                        else
                            conv[stnpair] = [cor(ref[:,1], stackData[stnpair][:,1])]
                        end
                    end
                elseif metric=="snr"
                    # compute peak signal-to-noise ratio of increasingly long stacks
                    psnr = maximum(snr(stackData[stnpair], reference.fs))
                    if haskey(conv, stnpair)
                        push!(conv[stnpair], psnr)
                    else
                        conv[stnpair] = [psnr]
                    end
                end
            end
        end
        close(f_cur) # xcorrs for each time step
    end
    close(ref_f) # reference xcorrs

    # write convergence vectors to disk
    f_conv = jldopen(foname, "a+")
    for key in keys(conv)
        f_conv["$metric/$key"] = conv[key]
    end
    close(f_conv) # convergence file
end

function daily_stack(data::CorrData; threshold::Float64=0.30, iterations::Int64=3)
    # preliminary stack of all constituent cross-correlations
    stackedData = stack(data, allstack=true)
    # snr array (len(cross-correlation), numIterations)
    snrArr = Array{Float32, 2}(undef, size(data.corr)[1], iterations+1)
    snrArr[:, 1] = snr(data,smooth=false)

    removedList = Array{Int64, 1}(undef, iterations)

    # number of iterations to remove incoherent cross-correlations
    p=PlotlyJS.plot([NaN])
    for i=1:iterations
        trace = PlotlyJS.scatter(;y=stackedData.corr[:,1])
        addtraces!(p, trace)

        # array of correlation coefficients for each constituent cross-correlation window
        ccList = Array{Float32,1}(undef,size(data.corr)[2])
        # iterate over each windowed cross-correlation
        for j=1:size(data.corr)[2]
            ccList[j] = cor(stackedData.corr[:, 1], data.corr[:, j])
        end

        #t=PlotlyJS.scatter(;y=ccList, mode="markers")
        #p1=PlotlyJS.plot([t])
        #display(p1)
        #readline()
        # find cross-correlations that fall below the cc threshold
        # on the next iteration, do not include these in the reference stack
        good_fit = findall(x->(x>=threshold), ccList)

        removedList[i] = size(data.corr)[2] - length(good_fit)

        # get SNR for cross-correlations that exceed the correlation coefficient threshold
        tempData = copy(data)
        tempData.corr = tempData.corr[:, good_fit]

        snrArr[:, i+1] = snr(tempData, smooth=false)

        stackedData = stack(tempData, allstack=true)
    end
    trace=PlotlyJS.scatter(;y=stackedData.corr[:,1])
    addtraces!(p, trace)
    deletetraces!(p, 1)
    display(p)
    readline()
    return stackedData, snrArr, removedList
end

f=jldopen("/Users/jared/SCECintern2019/data/xcorrs/BPnetwork_2003/BPnetwork_2003_xcorrs.2003.1.T00:00:00.jld2", "r")
data=f["2003.1.T00:00:00/BP.CCRB..BP1.BP.GHIB..BP1"]
stackedData, snrArr, removedList = daily_stack(data, threshold=0.3)
p=PlotlyJS.plot(snrArr)
display(p)
readline()
close(f)

p=PlotlyJS.plot(removedList)
display(p)
readline()
