include("partition.jl")
include("io.jl")

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
function convergence(basefiname::String, refname::String, foname::String; ntimes::Int64=-1, metric::String="cc", skipoverlap::Bool=false, slice::Union{Bool, Float64, Array{Float64,1}}=false, thresh::Float64=-1)
    # input file contains metadata (stationlist, timestamplist)
    f = jldopen(basefiname*".jld2")
    tslist = f["info/timestamplist"] # time stamps
    close(f) # input file

    # read reference xcorrs
    f_ref = jldopen(refname)

    # dictionary to keep track of convergence vectors for each station pair
    conv = Dict()
    # dictionary to keep track of stacked data
    stackData = Dict{String, Array{Float64,2}}()

    if ntimes==-1
        ntimes=length(tslist)
    end
    nWins = Dict()
    nGood = Dict()
    # iterate over timestamps
    for tstamp in tslist[1:ntimes]
        println("Processing $tstamp")
        # read xcorrs for each time step
        f_cur = jldopen(basefiname*".$tstamp.jld2")
        stnpairs = sort(keys(f_cur[tstamp]))
        stnpairs = ["BP.CCRB..BP1.BP.SMNB..BP1"]

        # iterate over station pairs
        for stnpair in stnpairs
            # load unstacked cross-correlations
            data = try f_cur["$tstamp/$stnpair"] catch; continue end

            # load reference cross-correlation for the current station pair
            reference = f_ref[stnpair]

            # optionally slice reference to keep only coda or only ballistic
            ref = corr_slice(reference, 1:size(reference.corr, 2), slice)

            if haskey(nWins, stnpair)
                push!(nWins[stnpair], size(data.corr, 2))
            else
                nWins[stnpair] = [size(data.corr, 2)]
            end

            datacopy, cList = selective_stacking(data, reference, threshold=thresh)
            good_inds = findall(x->(x>thresh), cList)
            data.corr = data.corr[:, good_inds]

            if haskey(nGood, stnpair)
                push!(nGood[stnpair], length(good_inds))
            else
                nGood[stnpair] = [length(good_inds)]
            end

            # iterate over windowed cross-correlations
            for i=1:size(data.corr)[2]
                if skipoverlap i=2*i-1 end # assumes cc_step is 1/2 cc_len

                # slice current cross-correlation if necessary
                data_cur = corr_slice(data, i:i, slice)

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
                    if haskey(conv, stnpair)
                        # compute psnr for both positive and negative coda
                        if typeof(slice)==Array{Float64,1}
                            neg_coda = findall(x->(x%2!=0), collect(1:size(stackData[stnpair], 2)))
                            pos_coda = findall(x->(x%2==0), collect(1:size(stackData[stnpair], 2)))
                            neg_coda_snr = maximum(snr(stackData[stnpair][:, neg_coda], reference.fs))
                            pos_coda_snr = maximum(snr(stackData[stnpair][:, pos_coda], reference.fs))

                            push!(conv[stnpair], [neg_coda_snr, pos_coda_snr])
                        else
                            push!(conv[stnpair], maximum(snr(stackData[stnpair], reference.fs)))
                        end
                    else
                        if typeof(slice)==Array{Float64,1}
                            neg_coda = findall(x->(x%2!=0), collect(1:size(stackData[stnpair], 2)))
                            pos_coda = findall(x->(x%2==0), collect(1:size(stackData[stnpair], 2)))
                            neg_coda_snr = maximum(snr(stackData[stnpair][:, neg_coda], reference.fs))
                            pos_coda_snr = maximum(snr(stackData[stnpair][:, pos_coda], reference.fs))

                            conv[stnpair] = [[neg_coda_snr, pos_coda_snr]]
                        else
                            conv[stnpair] = [maximum(snr(stackData[stnpair], reference.fs))]
                        end
                    end
                end
            end
        end
        close(f_cur) # xcorrs for each time step
    end
    close(f_ref) # reference xcorrs

    # write convergence vectors to disk
    f_conv = jldopen(foname, "a+")
    for key in keys(conv)
        f_conv["$metric/$key"] = conv[key]
        f_conv["$metric/nWin/$key"] = nWins[key]
        f_conv["$metric/nGood/$key"] = nGood[key]
    end

    close(f_conv) # convergence file
end

function corr_slice(data::CorrData, range::UnitRange{Int64}, slice::Union{Bool, Float64, Array{Float64,1}}=false)
    if typeof(slice)==Float64
        # time vector
        tvec = -data.maxlag:1/data.fs:data.maxlag
        # find all times within [-slice,slice] time
        t_inds = findall(x->(x>=-slice && x<=slice), tvec)
        # slice reference (maintain 2d array)
        d = data.corr[t_inds, range]
    elseif typeof(slice)==Array{Float64,1}
        # convert startlag/endlag[s] to startlag/windowlength[samples]
        win_len = Int(diff(slice)[1] * data.fs)
        startlag = Int(slice[1] * data.fs)
        # partition reference to keep only the coda
        d = partition(data.corr[:, range], startlag, win_len)
    else
        println("Proceeding without slicing.")
        d = data.corr
    end
    return d
end
