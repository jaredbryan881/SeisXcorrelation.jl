using SeisIO, SeisNoise, JLD2, PlotlyJS, StatsBase, Sockets, ORCA, Statistics

include("../../src/pairing.jl")
include("../../src/reference.jl")

function plot_seismograms(finame::String, stn::String; norm_factor=nothing, sparse::Int64=1, foname::String="", show::Bool=true, timeslice::Array{Int64,1}=[1, -1])
    f = jldopen(finame)
    timestamplist = f["info/DLtimestamplist"]

    if timeslice[2]==-1
        timeslice = [1, length(timestamplist)]
    end

    # set up plot
    layout = Layout(width=1200, height=2000,
                    xaxis_title="Time [s]",
                    yaxis_title="Days since $(timestamplist[1])",
                    font=attr(size=18),
                    showlegend=false,
                    title="$stn")
    p = PlotlyJS.plot([NaN], layout)

    for (titer, time) in enumerate(timestamplist[timeslice[1]:timeslice[2]])
        ts = f["$time/$stn"]
        if ts.misc["dlerror"] == 1 continue end
        # extract SeisChannel if necessary
        if typeof(ts)==SeisData ts=ts[1] end

        taxis = collect(0:length(ts.x)-1) ./ ts.fs

        if norm_factor == nothing
            norm = maximum(ts.x)
        else
            norm = norm_factor
        end

        trace = scatter(;x=taxis[1:sparse:end], y=((ts.x ./ (norm)) .+ titer .- 1)[1:sparse:end], mode="lines", line_color="rgb(0,0,0)", name=time)
        addtraces!(p, trace)
    end
    deletetraces!(p, 1)
    if foname !== ""
        savefig(p, foname)
    end
    if show
        display(p)
        readline()
    end

    close(f)
end

function plot_seismograms(finame::String; basefoname::String="", show::Bool=false, sparse::Int64=1)
    f=jldopen(finame)
    stlist = f["info/stationlist"]
    for stn in stlist
        println("Plotting $stn")
        plot_seismograms(finame, stn, foname="$basefoname.$stn.pdf", show=show, sparse=sparse)
    end
end

function plot_xcorrs(basefiname::String, stn1::String, stn2::String; reference::Union{Bool, String}=false, type::String="wiggles", foname::String="", show::Bool=true, filter::Union{Bool, Array{Float64,1}}=false, phase_smoothing::Float64=0., timeslice::Array{Int64,1}=[1,-1], stack::String="selective")
    f = jldopen(basefiname*".jld2")
    timestamplist = f["info/timestamplist"]

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    yaxis_title="Days since $(timestamplist[1])",
                    font=attr(size=18),
                    showlegend=false,
                    title="$stn1.$stn2")
    p = PlotlyJS.plot([NaN], layout)

    stnpair = stn1*"."*stn2
    stnpairrev = stn2*"."*stn1

    if type=="heatmap" xcorr_heat = []; tscounter=0 end

    if timeslice[2]==-1
        timeslice=[1,length(timestamplist)]
    end

    for (titer, time) in enumerate(timestamplist[timeslice[1]:timeslice[2]])
        f_cur = jldopen(basefiname*".$time.jld2")

        # load cross-correlation or its reverse
        if stnpair ∈ keys(f_cur[time])
            xcorr = f_cur["$time/$stnpair"]
            if stack=="selective"
                if reference!=false
                    f_ref=jldopen(reference)
                    ref = f_ref[stnpair]
                    xcorr, nRem = selective_stacking(xcorr, ref, threshold=0.3)
                    close(f_ref)
                else
                    xcorr, nRem = selective_stacking(xcorr, threshold=0.3)
                end
            else
                stack!(xcorr, allstack=true, phase_smoothing=phase_smoothing)
            end

        elseif stnpairrev ∈ keys(f_cur[time])
            xcorr = f_cur["$time/$stnpairrev"]
            if stack=="selective"
                if reference!=false
                    f_ref=jldopen(reference)
                    ref = f_ref[stnpairrev]
                    xcorr, nRem = selective_stacking(xcorr, ref, threshold=0.3)
                    close(f_ref)
                else
                    xcorr, nRem = selective_stacking(xcorr, threshold=0.3)
                end
            else
                stack!(xcorr, allstack=true, phase_smoothing=phase_smoothing)
            end
            xcorr.corr = reverse(xcorr.corr, dims=1)

        else
            xcorr = CorrData()
            xcorr.corr = zeros(4001,1)
        end

        if filter!=false
            xcorr.corr[:,1] = bandpass(xcorr.corr[:,1], filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
        end

        # make lags visible outside of loop
        global lags = -xcorr.maxlag:1/xcorr.fs:xcorr.maxlag

        norm_factor = maximum(abs.(xcorr.corr[:]))

        try
            if type=="wiggles"
                # add trace for wiggle plot
                trace = scatter(;x=lags, y=xcorr.corr[:, 1] ./ (2 * norm_factor) .+ titer .- 1, mode="lines", line_color="black", name=time)
                addtraces!(p, trace)

            elseif type=="heatmap"
                # append to array for heat map. We cannot plot before it is fully filled
                append!(xcorr_heat, xcorr.corr[:, 1]./ (norm_factor))
                tscounter+=1
            end
        catch y
            println(y)
            println("Plotting failed")
        end
        close(f_cur)
    end

    if type=="heatmap"
        # dont plot if no cross-correlations were found. This would give div by 0
        if tscounter==0 println("No cross-correlations to plot. Exiting."); exit() end
        xcorr_heat = reshape(xcorr_heat, convert(Int64, length(xcorr_heat)/tscounter), tscounter)
        trace = PlotlyJS.heatmap(x=lags, z=xcorr_heat)
        addtraces!(p, trace)
    end
    deletetraces!(p, 1)

    if foname !== ""
        savefig(p, foname)
    end
    if show
        display(p)
        readline()
    end

    close(f)
end

function plot_xcorrs(basefiname::String, basefoname::String; show::Bool=false, type::String="heatmap")
    f=jldopen(basefiname*".jld2")
    stlist = f["info/stationlist"]
    stniter=1
    for stn1 in stlist
        for stn2 in stlist[stniter:end]
            println("Processing $stn1--$stn2")
            plot_xcorrs(basefiname, stn1, stn2; type=type, foname=basefoname*".$stn1.$stn2.png", show=show)
        end
        stniter+=1
    end
end

function plot_reference(finame::String, stn1::String, stn2::String)
    f = jldopen(finame)
    stname = stn1*"."*stn2
    data = f[stname]

    lags = -data.maxlag:1/data.fs:data.maxlag

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    font =attr(size=18),
                    showlegend=true,
                    title="$stname",
                    paper_bgcolor="rgba(0,0,0,0)",
                    plot_bgcolor="rgba(0,0,0,0)")

    p = PlotlyJS.plot([NaN], layout)
    trace1 = scatter(;x=lags, y=data.corr[:, 1], mode="lines", line_color="rgb(0,0,0)")
    addtraces!(p, trace1)
    deletetraces!(p, 1)
    display(p)
    readline()
end

function plot_convergence(finame::String, stn1::String, stn2::String; maxlag::Float64=100.0, fs::Float64=20.0)
    f = jldopen(finame)
    stname = stn1*"."*stn2
    data = f[stname]

    lags = -maxlag:1/fs:maxlag

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="# of daily stacked cross-correlations",
                    font=attr(size=18),
                    showlegend=true,
                    title="$stname RMS convergence")

    p = PlotlyJS.plot([NaN], layout)

    trace1 = scatter(;x=collect(1:length(data))./46, y=data, mode="lines", line_color="rgb(0,0,0)")
    addtraces!(p, trace1)
    deletetraces!(p, 1)
    display(p)
    readline()
end

function plot_convergence(finame::String; sparse::Int64=1)
    f=jldopen(finame)
    lags = -100.0:0.05:100.0

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="# of daily stacked cross-correlations",
                    font=attr(size=18),
                    showlegend=true,
                    title="RMS Convergence")
    p = PlotlyJS.plot([NaN], layout)
    for stpair in keys(f)[1:200]
        (net1, stn1, comp1, net2, stn2, comp2) = filter!(x->x≠"", split(stpair, "."))
        if comp1==comp2 && stn1!=stn2 && comp1=="BP3"
            color="rgba(150,0,0,1.0)"
        elseif stn1==stn2
            continue
        else
            color="rgba(0,0,0,0.1)"
        end
        data = f[stpair]
        trace = scatter(;x=(collect(1:length(data))[1:sparse:end])./46, y=data[1:sparse:end], mode="lines", line_color=color, name=stpair)
        addtraces!(p, trace)
    end
    deletetraces!(p, 1)
    display(p)
    readline()
end

function count_zeros(ref::Union{Array{Float32,1}, Array{Float64,1}}, data::Union{Array{Float32,1}, Array{Float64,1}}; percent::Bool=true)
    # get number of zeros in two given arrays, and compute the difference.
    # in the expected use case, this finds the amount of data zeroed by earthquake removal
    nZeros = count(x->x==0, data)
    nZerosRef = count(x->x==0, ref)
    rem = nZeros-nZerosRef
    # convert number of removed zeros to percent of total data length
    if percent
        rem /= length(ref)
        rem *= 100
    end
    return rem
end

function plot_removal(orig_finame::String, finame::String; kurtlist::Array{Float64,1}, staltalist::Array{Float64,1})
    # open original input file
    f_orig = jldopen(orig_finame*".jld2")
    # time stamps and stations
    tslist = f_orig["info/DLtimestamplist"]
    stlist = f_orig["info/stationlist"]

    # initialize 2D array to store number of removed data points
    removed = zeros(length(staltalist), length(kurtlist), length(stlist))
    # initialize 3D array to store peak amplitude
    amps = zeros(length(staltalist), length(kurtlist), length(tslist), length(stlist))
    # initialize 3D array to store average amplitude
    avgs = zeros(length(staltalist), length(kurtlist), length(tslist), length(stlist))

    for (i, stalta) in enumerate(staltalist)
        for (j, kurt) in enumerate(kurtlist)
            println("Processing stalta: $stalta, kurtosis: $kurt")
            # open file for each stalta threshold and kurtosis threshold
            f=jldopen(finame*"_$(stalta)_$(kurt).jld2")
            # count number of zeros for all time stamps
            for (t,ts) in enumerate(tslist)
                for (s, stn) in enumerate(stlist)
                    orig_data = f_orig["$ts/$stn"].x
                    data = f["$ts/$stn"].x
                    removed[i,j,s] += count_zeros(orig_data, data)
                    amps[i,j,t,s] = maximum(abs.(data))
                    avgs[i,j,t,s] = mean(data)
                end
            end
            close(f)

            # average percent removal over days
            removed[i,j,:] ./= length(tslist)
        end
    end
    close(f_orig)

    for (s, stn) in enumerate(stlist)
        # set up plot
        layout = Layout(width=1200, height=800,
                        xaxis_title="STA/LTA Threshold",
                        yaxis_title="Kurtosis Threshold",
                        font=attr(size=18),
                        showlegend=true,
                        title="$stn % removal")

        p = PlotlyJS.plot([NaN], layout)
        trace = PlotlyJS.heatmap(x=staltalist, y=kurtlist, z=removed[:,:,s])
        addtraces!(p, trace)
        deletetraces!(p, 1)
        display(p)
        readline()
    end
end

function plot_dvv(finame::String, stn1::String, stn2::String)
    f=jldopen(finame)
    dvv = f["$stn1.$stn2"]
    tslist = f["info/timestamplist"]

    layout = Layout(width=1200, height=800,
                    xaxis_title="Days since $(tslist[1])",
                    yaxis_title="dv/v (%)",
                    font=attr(size=18),
                    showlegend=false,
                    title="$stn1--$stn2")

    p = PlotlyJS.plot([NaN], layout)
    trace = PlotlyJS.scatter(;x=collect(1:length(tslist)), y=dvv, mode="lines+markers", line_color="rgba(0,0,0,1)")
    addtraces!(p, trace)
    deletetraces!(p, 1)
    display(p)
    readline()
end

function plot_dvv(finame::String, stack=false)
    f=jldopen(finame)
    tslist = f["info/timestamplist"]

    layout = Layout(width=1200, height=800,
                    xaxis_title="Days since $(tslist[1])",
                    yaxis_title="dv/v (%)",
                    font=attr(size=18),
                    showlegend=false)
    p = PlotlyJS.plot([NaN], layout)

    if stack
        dvvfull = zeros(31)
        count = 0
    end

    for stpair in keys(f)[1:780]
        (net1, stn1, comp1, net2, stn2, comp2) = filter!(x->x≠"", split(stpair, "."))
        if stpair == "BP.CCRB..BP1.BP.SMNB..BP1"
            color="rgba(150,0,0,1)"
        elseif stn1==stn2
            continue
        elseif stn1!=stn2 && comp1!=comp2
            continue
        else
            color = "rgba(0,0,0,0.2)"
        end

        if stack
            dvv = f[stpair]
            dvvfull .+= dvv
            count += 1
        else
            trace = PlotlyJS.scatter(;x=collect(1:length(tslist)), y=dvv, mode="lines+markers", line_color=color)
            addtraces!(p,trace)
        end
    end

    if stack
        trace = PlotlyJS.scatter(;x=collect(1:length(tslist)), y=dvvfull./count, mode="lines+markers", line_color="rgba(0,0,0,1)")
        addtraces!(p, trace)
    end

    close(f)
    deletetraces!(p,1)
    display(p)
    readline()
end

function plot_dvv(basefiname::String, methods::Array{String,1}, percent::Array{Bool,1})
    f1=jldopen(basefiname*"_$(methods[1]).jld2")
    tslist = f1["info/timestamplist"]
    close(f1)
    layout = Layout(width=1200, height=800,
                    xaxis_title="Days since $(tslist[1])",
                    yaxis_title="dv/v (%)",
                    yaxis_range=[-0.1, 0.1],
                    font=attr(size=18),
                    showlegend=true)
    p = PlotlyJS.plot([NaN], layout)

    for (i,method) in enumerate(methods)
        f=jldopen(basefiname*"_$(methods[i]).jld2")

        dvvfull = zeros(31)
        count=0

        for stpair in keys(f)[1:780]
            (net1, stn1, comp1, net2, stn2, comp2) = filter!(x->x≠"", split(stpair, "."))
            if stpair == "BP.CCRB..BP1.BP.SMNB..BP1"
                color="rgba(150,0,0,1)"
            elseif stn1==stn2
                continue
            elseif stn1!=stn2 && comp1!=comp2
                continue
            else
                color = "rgba(0,0,0,0.2)"
            end

            dvv = f[stpair]
            dvvfull .+= dvv
            count += 1
        end

        if percent[i]!=true dvvfull .*= 100 end
        trace = PlotlyJS.scatter(;x=collect(1:length(tslist)), y=dvvfull./count, mode="lines+markers", line_color="rgba(0,0,0,1)", name="$method")
        addtraces!(p, trace)
        close(f)
    end
    deletetraces!(p, 1)
    display(p)
    readline()
end
