using SeisIO, Noise, JLD2, PlotlyJS

# This is a script to test the functionality of the MWCS functionality of Noise.jl
# This assumes synthetic data generated via makeStretchData.jl

finame = "StretchedDataSmall"
valData = jldopen(finame*".jld2")
u1 = valData["u0"]
u2 = valData["u1"]
st = valData["st"]
dt = valData["dt"]

# MWCS input parameters
InputDict = Dict( "freqmin"            => 0.001,
                  "freqmax"            => 9.9,
                  "fs"                 => 20.0,
                  "mintimelag"         => 0.0,
                  "dtt_width"          => 30.0,
                  "window_length"      => 10.0,
                  "window_step"        => 10.0,
                  "dtt_lag"            => "static",
                  "dtt_v"              => 1.0,
                  "dtt_minlag"         => 0.0,
                  "dtt_sides"          => "Both",
                  "smoothing_half_win" => 3)
dist = 0.0
true_dtt = 0.0001
# compute dt
time_axis, dt, err, coh = mwcs(u1,
                               u2,
                               InputDict["freqmin"],
                               InputDict["freqmax"],
                               InputDict["fs"],
                               InputDict["mintimelag"],
                               InputDict["window_length"],
                               InputDict["window_step"],
                               InputDict["smoothing_half_win"])
# compute WLS regression for dt/t
m, em, a, ea, m0, em0 = mwcs_dvv(time_axis,
                                 dt,
                                 err,
                                 coh,
                                 InputDict["dtt_lag"],
                                 dist/1000,
                                 InputDict["dtt_v"],
                                 InputDict["dtt_minlag"],
                                 InputDict["dtt_width"],
                                 InputDict["dtt_sides"])

# Plot dt vs t, error bars, true dt vs t, and the best fit line
trace1 = scatter(;x=time_axis, y=dt, mode="lines+markers", marker_color=coh, line_color="black", name="Calculated dt/t: $((round(m, sigdigits=3))*100)%")
true_shift = scatter(;x=time_axis, y=true_dtt*time_axis, mode="lines", line_color="blue", name="True dt/t: $(true_dtt*100)%")
trace2 = scatter(;x=time_axis, y=(m*time_axis) .+ a, mode="lines", line_color="green", name="Linear Regression")
plots = [trace1, true_shift, trace2]
# plot error bars from mwcs()
err1 = dt .+ err
err2 = dt .- err
for i=1:length(dt)
    t = scatter(;x=[time_axis[i], time_axis[i]], y=[err1[i], err2[i]], mode="lines", line_color="black", name=nothing)
    push!(plots, t)
end
layout = Layout(;title="MWCS dt/t", xaxis_title="t (s)", yaxis_title="dt (s)")
p = PlotlyJS.plot(plots, layout)
display(p)
readline()
close(valData)

# Report performance
println("True: $(-true_dtt*100)%")
println("Estimated: $(-m*100)%")
error = -(m - true_dtt)*100
println("Error: $error")
