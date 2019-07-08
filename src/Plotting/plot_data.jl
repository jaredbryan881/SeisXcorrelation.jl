using SeisIO, SeisNoise, JLD2, PlotlyJS, StatsBase, Sockets, ORCA

function plot_seismograms(finame::String, stn::String; sparse::Int64=1, foname::String="", show::Bool=true)
    f = jldopen(finame)
    timestamplist = f["info/DLtimestamplist"]

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Time [s]",
                    yaxis_title="Days",
                    font=attr(size=18),
                    showlegend=true,
                    title="$stn")
    p = PlotlyJS.plot([NaN], layout)

    for (titer, time) in enumerate(timestamplist)
        ts = f["$time/$stn"]
        if typeof(ts)==SeisData ts=ts[1] end

        taxis = collect(0:length(ts.x)-1) ./ ts.fs
        norm_factor = maximum(ts.x)
        trace = scatter(;x=taxis[1:sparse:end], y=((ts.x ./ norm_factor) .+ titer .- 1)[1:sparse:end], mode="lines", linecolor="rgb(0,0,0)", name="$time")
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

function plot_xcorrs(basefiname::String, stn1::String, stn2::String; type::String="wiggles", maxlag::Float64=100.0, dt::Float64=0.05, foname::String="", show::Bool=true)
    f = jldopen(basefiname*".jld2")
    timestamplist = f["info/timestamplist"]

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    font=attr(size=18),
                    showlegend=true,
                    title="$stn1.$stn2")
    p = PlotlyJS.plot([NaN], layout)

    stnpair = stn1*"."*stn2
    stnpairrev = stn2*"."*stn1

    if type=="heatmap" xcorr_heat = [] end

    for (titer, time) in enumerate(timestamplist)
        f_cur = jldopen(basefiname*".$time.jld2")
        xcorr = try
            xcorr = f_cur["$time/$stnpair"]
            stack!(xcorr, allstack=true)
            xcorr
        catch y;
            xcorr = f_cur["$time/$stnpairrev"]
            stack!(xcorr, allstack=true)
            xcorr.corr = reverse(xcorr.corr, dims=1)
            xcorr
        end
        lags = -xcorr.maxlag:1/xcorr.fs:xcorr.maxlag

        norm_factor = maximum(xcorr.corr[:])

        try
            if type=="wiggles"
                trace = scatter(;x=lags, y=xcorr.corr[:,1] ./ (2 * norm_factor) .+ titer .- 1, mode="lines", linecolor="rgb(0,0,0)", name="$time")
                addtraces!(p, trace)
            elseif type=="heatmap"
                append!(xcorr_heat, xcorr.corr[:, 1]./ (norm_factor))
            end
        catch y
            println(y)
            println("Plotting failed")
            close(f_cur)
        end
        close(f_cur)
    end

    if type=="heatmap"
        lags = -maxlag:dt:maxlag
        xcorr_heat = reshape(xcorr_heat, convert(Int64, length(xcorr_heat)/length(timestamplist)), length(timestamplist))
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
    trace1 = scatter(;x=lags, y=data[:, 1], mode="lines", linecolor="rgb(0,0,0)")
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

    trace1 = scatter(;x=1:length(data), y=data, mode="lines", linecolor="rgb(0,0,0)")
    addtraces!(p, trace1)
    deletetraces!(p, 1)
    display(p)
    readline()
end
