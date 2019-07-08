using PlotlyJS, JLD2, ORCA

function plot_error(finame::String, type::String; show::Bool=true, foname::String="")
    f = jldopen(finame)
    dvVlist = f["info/dvV"]
    dvVlist = dvVlist[1:5:end]
    noiselist = f["info/noise"]

    errors = zeros(length(dvVlist), length(noiselist))
    # assemble matrix for heatmap
    for (i,dvV) in enumerate(dvVlist)
        for (j,noise) in enumerate(noiselist)
            errors[i,j] = f["$type/$dvV.$noise"]
        end
    end

    # set up plot
    layout = Layout(width=1200, height=800,
                    xaxis_title="dv/v",
                    yaxis_title="Noise level",
                    font=attr(size=18),
                    showlegend=true,
                    title="Stretching Error")
    p = PlotlyJS.plot([NaN], layout)
    trace = PlotlyJS.heatmap(x=dvVlist .* 100, y=noiselist .* 100, z=errors)
    addtraces!(p, trace)
    deletetraces!(p, 1)

    if show==true
        display(p)
        readline()
    end
    if foname !== ""
        savefig(p, foname)
    end
end

plot_error("verificationDataError_DTW.jld2", "realData", foname="DTWError_stretchXcorr.pdf")
