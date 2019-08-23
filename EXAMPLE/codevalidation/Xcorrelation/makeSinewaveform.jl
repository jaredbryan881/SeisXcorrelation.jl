using PlotlyJS, FileIO, JLD2, SeisIO

# This script makes two time_shifted or stretched traces for validation of SeisXcorrelation

#----------------------------------#
dt     = 0.01; # sampling time of waveform
# control time shift
tarrival_1 = 30 #wave arrival time at station 1
tarrival_2 = 30.05 #wave arrival time at station 2
# control relative frequency
f1 = 0.1 # source sine frequency Hz
f2 = f1*1.01

npts = 6000; # length of time series

foname = "./sinStretchTest.jld2"
#----------------------------------#
#Period of sine
T1 = 1/f1
T2 = 1/f2

npts_sine1 = round(Int64, T1/dt)
npts_sine2 = round(Int64, T2/dt)

# make waveform time series
s1   = sin.(LinRange(0, 2*pi, npts_sine1));
s2   = sin.(LinRange(0, 2*pi, npts_sine2));

u1 = zeros(Float64, npts)
u2 = zeros(Float64, npts)

tp1 = round(Int64, tarrival_1/dt)
tp2 = round(Int64, tarrival_2/dt)

u1[tp1:tp1+npts_sine1-1] = s1
u2[tp2:tp2+npts_sine2-1] = s2

tvec  = collect(( 0:npts - 1 ) .* dt);

Δtlag_true = tarrival_2 - tarrival_1 # true arrival time, which should be equal to the time at maximam coherency of cross-correlation function

# package generated signal in SeisChannel and save
t=[1 0; npts 0]
Ch1 = SeisChannel(name="BP.u1u1.fst", id="BP.u1u1.fst", x=u1, fs=1/dt, t=t)
Ch1.misc["dlerror"] = 0
Ch2 = SeisChannel(name="BP.u2u2.snd", id="BP.u2u2.snd", x=u2, fs=1/dt, t=t)
Ch2.misc["dlerror"] = 0
jldopen(foname, "a+") do file
    file["info/DLtimestamplist"] = ["time1"]
    file["info/stationlist"] = ["BP.u1u1.fst", "BP.u2u2.snd"]
    file["info/DL_time_unit"] = npts*dt
    file["info/f1"] = f1
    file["info/f2"] = f2
    file["info/Δtlag_true"] = Δtlag_true
    file["time1/BP.u1u1.fst"] = Ch1
    file["time1/BP.u2u2.snd"] = Ch2
end

#plot curve
function lineplot1()
    tr1 = scatter(;x=tvec, y=u1, mode="lines+markers", name="u1")
    tr2 = scatter(;x=tvec, y=u2, mode="lines+markers", name="u2")
    layout = Layout(
        xaxis=attr(title="u1"),
        yaxis=attr(title="u2"),
        title="True time shift = $Δtlag_true [s]"
    )
    p = plot([tr1, tr2], layout)
    display(p)
    readline()
end
