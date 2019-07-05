using FileIO, SeisIO, SeisyNoise, JLD2, PlotlyJS

function plot_data(finame, stname)
    f = jldopen(finame)
    timestamplist = f["info/timestamplist"][1:end-1]

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    font =attr(size=18),
                    showlegend=true,
                    title="$stname")
    p = PlotlyJS.plot([NaN], layout)

    day = 0
    for t in timestamplist[1:end-1]
        key = "$t/$stname"
        ts = try f[key] catch; continue end
        lags = -ts.maxlag:1/ts.fs:ts.maxlag

        stack!(ts, allstack=true)

        norm_factor = maximum(ts.corr[:])

        try
            trace1 = scatter(;x=lags, y=xcorr[:,1] ./ (2 * norm_factor) .+ day, mode="lines", linecolor="rgb(0,0,0)", name="$t")
            addtraces!(p, trace1)
        catch
            println("Plotting failed")
            close(f)
            exit()
        end
        day = day + 1
    end
    deletetraces!(p, 1)
    display(p)
    readline()

    close(f)
end

function plot_reference()
    finame_weq="../reference_xcorr.jld2"
    finame_neq="../reference_xcorr_neq.jld2"

    fweq = jldopen(finame_weq)
    fneq = jldopen(finame_neq)

    stname = "BP.CCRB..BP1.BP.EADB..BP1"

    data_weq = fweq[stname]
    data_neq = fneq[stname]

    lags = -100.0:1/20.0:100.0

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    font =attr(size=18),
                    showlegend=true,
                    title="$stname")

    p = PlotlyJS.plot([NaN], layout)

    trace1 = scatter(;x=lags, y=data_weq[:, 1], mode="lines", linecolor="rgb(0,0,0)", name="Raw Data")
    #addtraces!(p, trace1)
    trace2 = scatter(;x=lags, y=data_neq[:, 1], mode="lines", linecolor="rgb(0,0,0)", name="Earthquakes Removed")
    addtraces!(p, trace2)
    deletetraces!(p, 1)
    display(p)
    readline()
end

function plot_convergence()
    finame_weq="../rms_weq_wol.jld2"
    finame_neq="../rms_neq_wol.jld2"

    fweq = jldopen(finame_weq)
    fneq = jldopen(finame_neq)

    stname = "BP.CCRB..BP1.BP.EADB..BP1"

    lags = -100.0:1/20.0:100.0
    data_weq = fweq[stname]
    data_neq = fneq[stname]

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Stacked cross-correlations",
                    font =attr(size=18),
                    showlegend=true,
                    title="$stname")

    p = PlotlyJS.plot([NaN], layout)

    trace1 = scatter(;x=1:length(data_weq), y=data_weq, mode="lines", linecolor="rgb(0,0,0)", name="Raw Data")
    addtraces!(p, trace1)
    trace2 = scatter(;x=1:length(data_weq), y=data_neq, mode="lines", linecolor="rgb(0,0,0)", name="Earthquakes Removed")
    #addtraces!(p, trace2)
    deletetraces!(p, 1)
    display(p)
    readline()
end
