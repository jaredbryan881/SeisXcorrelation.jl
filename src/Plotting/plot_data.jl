using SeisIO, SeisNoise, JLD2, PlotlyJS, StatsBase, Sockets, ORCA

function plot_seismograms(finame::String, stn::String; norm_factor=nothing, sparse::Int64=1, foname::String="", show::Bool=true)
    f = jldopen(finame)
    timestamplist = f["info/DLtimestamplist"]
    timestamplist=timestamplist[1:end-1]

    # set up plot
    layout = Layout(width=1200, height=1500,
                    xaxis_title="Time [s]",
                    yaxis_title="Days since $(timestamplist[1])",
                    font=attr(size=18),
                    showlegend=false,
                    title="$stn")
    p = PlotlyJS.plot([NaN], layout)

    for (titer, time) in enumerate(timestamplist)
        ts = f["$time/$stn"]
        if typeof(ts)==SeisData ts=ts[1] end

        taxis = collect(0:length(ts.x)-1) ./ ts.fs

        if norm_factor == nothing
            norm_factor = maximum(ts.x)
        end

        trace = scatter(;x=taxis[1:sparse:end], y=((ts.x ./ norm_factor) .+ titer .- 1)[1:sparse:end], mode="lines", line_color="rgb(0,0,0)", name=time)
        addtraces!(p, trace)
    end
    deletetraces!(p, 1)
    if foname !== ""
        savefig(p, foname)
    end
    if show == true
        display(p)
    end
    readline()

    close(f)
end

function plot_xcorrs(basefiname::String, stn1::String, stn2::String; type::String="wiggles", foname::String="", show::Bool=true)
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

    for (titer, time) in enumerate(timestamplist)
        f_cur = jldopen(basefiname*".$time.jld2")

        if stnpair ∈ keys(f_cur[time])
            xcorr = f_cur["$time/$stnpair"]
            stack!(xcorr, allstack=true)
        elseif stnpairrev ∈ keys(f_cur[time])
            xcorr = f_cur["$time/$stnpairrev"]
            stack!(xcorr, allstack=true)
            xcorr.corr = reverse(xcorr.corr, dims=1)
        else
            continue
        end

        global lags = -xcorr.maxlag:1/xcorr.fs:xcorr.maxlag

        norm_factor = maximum(xcorr.corr[:])

        try
            if type=="wiggles"
                trace = scatter(;x=lags, y=xcorr.corr[:,1] ./ (2 * norm_factor) .+ titer .- 1, mode="lines", line_color="black", name=time)
                addtraces!(p, trace)
            elseif type=="heatmap"
                append!(xcorr_heat, xcorr.corr[:, 1]./ (norm_factor))
                tscounter+=1
            end
        catch y
            println(y)
            println("Plotting failed")
            close(f_cur)
        end
        close(f_cur)
    end

    if type=="heatmap"
        if tscounter==0 println("No cross-correlations to plot. Exiting."); exit() end
        xcorr_heat = reshape(xcorr_heat, convert(Int64, length(xcorr_heat)/tscounter), tscounter)
        trace = PlotlyJS.heatmap(x=lags, z=xcorr_heat)
        addtraces!(p, trace)
    end
    deletetraces!(p, 1)
    if foname !== ""
        savefig(p, foname)
    end
    if show == true
        display(p)
    end
    readline()

    close(f)
end

function plot_reference(finame::String, stn1::String, stn2::String)
    f = jldopen(finame)
    stname = stn1*"."*stn2
    data = f[stname]

    lags = -100.0:1/20.0:100.0

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    font =attr(size=18),
                    showlegend=true,
                    title="$stname")

    p = PlotlyJS.plot([NaN], layout)
    trace1 = scatter(;x=lags, y=data[:, 1], mode="lines", line_color="rgb(0,0,0)")
    addtraces!(p, trace1)
    deletetraces!(p, 1)
    display(p)
    readline()
end

function plot_convergence(finame::String, stn1::String, stn2::String)
    f = jldopen(finame)
    stname = stn1*"."*stn2
    data = f[stname]

    lags = -100.0:1/20.0:100.0

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="# of stacked cross-correlations",
                    font=attr(size=18),
                    showlegend=true,
                    title="$stname convergence")

    p = PlotlyJS.plot([NaN], layout)

    trace1 = scatter(;x=1:length(data), y=data, mode="lines", line_color="rgb(0,0,0)")
    addtraces!(p, trace1)
    deletetraces!(p, 1)
    display(p)
    readline()
end
