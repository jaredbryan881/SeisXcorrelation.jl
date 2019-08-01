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
                    close(f_cur)
                    close(ref_f)
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
