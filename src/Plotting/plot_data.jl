using FileIO, SeisIO, Noise, JLD2, PlotlyJS

function plot_data()
    finame="/Users/jared/SCECintern2019/SeisXcorrelation.jl/EXAMPLE/xcorr_BP/outputData/BPnetworkxcorr_weq.jld2"
    stnname="BP.EADB..BP1.BP.VCAB..BP1"
    f = jldopen(finame)
    timestamplist = f["info/timestamplist"][1:end-1]

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="Lags [s]",
                    font =attr(size=18),
                    showlegend=true,
                    title="$stnname")
    p = PlotlyJS.plot([NaN], layout)

    day = 0
    for t in timestamplist
        key = "$t/$stnname"
        ts = f[key]
        lags = -ts.maxlag:1/ts.fs:ts.maxlag

        stack!(ts, allstack=true)

        norm_factor = maximum(ts.corr[:])

        xcorr = bandpass(ts.corr[:, 1], 0.1, 0.9, ts.fs)
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

plot_data()
