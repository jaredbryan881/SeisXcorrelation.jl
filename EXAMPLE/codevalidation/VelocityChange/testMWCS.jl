using SeisIO, SeisNoise, JLD2, PlotlyJS

# This is a script to test the functionality of the MWCS functions of Noise.jl
# This assumes synthetic data generated via genSynth.jl

function testMWCS(finame::String, InputDict::Dict, windowlenlist::Array{Float64,1}, winsteplist::Array{Float64,1}, smoothlist::Array{Int64,1}, type::String)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    for dvV in dvVlist[1:10:end]
        for noiselvl in noiselist[1:1]
            for win_len in windowlenlist, win_step in windowsteplist, smooth in smoothlist
                signals = valData["$type/$dvV.$noiselvl"]
                u1 = signals[1]
                u2 = signals[2]

                dist = 0.0
                true_dtt = -dvV
                # compute dt
                time_axis, dt, err, coh = try
                    mwcs(u1,
                         u2,
                         InputDict["freqmin"],
                         InputDict["freqmax"],
                         InputDict["fs"],
                         InputDict["mintimelag"],
                         win_len,
                         win_step,
                         smooth)
                catch y
                    println("Failed on first step")
                    println(y)
                    continue
                end
                dt .*= 2π
                err .*= 2π
                # compute WLS regression for dt/t
                m, em, a, ea, m0, em0 = try
                    mwcs_dvv(time_axis,
                             dt,
                             err,
                             coh,
                             InputDict["dtt_lag"],
                             dist/1000,
                             InputDict["dtt_v"],
                             InputDict["dtt_minlag"],
                             InputDict["dtt_width"],
                             InputDict["dtt_sides"])
                catch y
                    println("Failed on second step")
                    println(y)
                    continue
                end

                # Report performance
                println("Smoothing: $smooth")
                println("Window_length: $win_len")
                println("True: $(-true_dtt*100)%")
                println("Estimated: $(-m0*100)%")
                error = abs((m0 - true_dtt)/true_dtt)*100
                println("Error: $error%")
                plot_dtt(time_axis, dt, coh, m0, 0, err, true_dtt)
            end
        end
    end
    close(valData)
end

function plot_dtt(time_axis, dt, coh, m, a, err, true_dtt)
    # Plot dt vs t, error bars, true dt vs t, and the best fit line
    trace1 = scatter(;x=time_axis, y=dt, mode="lines+markers", marker_color=coh, line_color="black", name="Calculated dt/t: $((round(m, sigdigits=3))*100)%")
    true_shift = scatter(;x=time_axis, y=true_dtt*time_axis, mode="lines", line_color="blue", name="True dt/t: $(true_dtt*100)%")
    trace2 = scatter(;x=time_axis, y=(m*time_axis) .+ a, mode="lines", line_color="green", name="Linear Regression")
    plots = [trace1, true_shift, trace2]

    # plot error bars from mwcs()
    err1 = dt .+ err
    err2 = dt .- err
    for i=1:length(dt)
        t = scatter(;x=[time_axis[i], time_axis[i]], y=[err1[i], err2[i]], mode="lines", line_color="black", showlegend=false)
        push!(plots, t)
    end

    layout = Layout(;title="MWCS dt/t", xaxis_title="t (s)", yaxis_title="dt (s)")
    p = PlotlyJS.plot(plots, layout)
    display(p)
    readline()
end

# MWCS input parameters
InputDict = Dict( "freqmin"            => 0.1,
                  "freqmax"            => 9.9,
                  "fs"                 => 20.0,
                  "mintimelag"         => 0.0,
                  "dtt_width"          => 200.0,
                  "window_length"      => 0,
                  "window_step"        => 0,
                  "dtt_lag"            => "static",
                  "dtt_v"              => 1.0,
                  "dtt_minlag"         => 0.0,
                  "dtt_sides"          => "both",
                  "smoothing_half_win" => 0)

finame = "verificationData.jld2"
windowlenlist = collect(40.:10.:60.)
windowsteplist = [10.0]
smoothlist = [7, 11, 15, 19]
type = "dampedSinusoid"
testMWCS(finame, InputDict, windowlenlist, windowsteplist, smoothlist, type)
