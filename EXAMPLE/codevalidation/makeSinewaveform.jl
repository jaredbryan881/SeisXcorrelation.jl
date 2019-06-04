using PlotlyJS, FileIO, JLD2

# This script makes two time_shifted traces for validation of SeisXcorrelation

#----------------------------------#
dt     = 0.05; # sampling time of waveform

tarrival_1 = 30 #wave arrival time at station 1
tarrival_2 = 40 #wave arrival time at station 2

f = 0.1 # source sine frequency Hz

npts = 1200; # length of time series
#----------------------------------#

T = 1/f #Period of sine

npts_sine = round(Int64, T/dt)

s0   = sin.(LinRange(0, 2*pi, npts_sine)); # make waveform time series

u1 = zeros(Float64, npts)
u2 = zeros(Float64, npts)

tp1 = round(Int64, tarrival_1/dt)
tp2 = round(Int64, tarrival_2/dt)

u1[tp1:tp1+npts_sine-1] = s0
u2[tp2:tp2+npts_sine-1] = s0

tvec  = collect(( 0:npts - 1 ) .* dt);

Δtlag_true = tarrival_2 - tarrival_1 # true arrival time, which should be equal to the time at maximam coherency of cross-correlation function

@save "./ValidationInput.jld2" dt u1 u2 Δtlag_true


#plot curve

function lineplot1()
    tr1 = scatter(;x=tvec, y=u1, mode="lines+markers", name="u1")
    tr2 = scatter(;x=tvec, y=u2, mode="lines+markers", name="u2")
    layout = Layout(
        xaxis=attr(title="u1"),
        yaxis=attr(title="u2"),
        title="True time shift = $Δtlag_true [s]"
    )
    plot([tr1, tr2], layout)
end

lineplot1()
