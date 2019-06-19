using PlotlyJS, Random, DSP, Dierckx, FileIO, JLD2

# This script makes the trace data to compare a step in the tine shifts

#ricker wavelet function
"""
ricker(omega)
returns Ricker wavelet

input:

    peak frequency   = 20 Hz
    sampling time    = 0.001 seconds
    number of points = 100;

output:
    s: ricker wavelet with length = number of points
    t: time
"""
function ricker(;f::Float64=20.0, n::Int64=100, dt::Float64=0.001, nargout::Bool=true)
    # Create the wavelet and shift in time if needed
    T = dt*(n-1);
    t = collect(0:dt:T);
    t0 = 1/f;
    tau = t.-t0;

    s = (1.0 .- tau.*tau.*f.^2*π.^2).*exp.(-tau.^2*π.^2*f.^2);

    if nargout
        trace = scatter(;x=t, y=s, mode="lines")
        layout = Layout(
            xaxis=attr(title="Time"),
            yaxis=attr(title="u"),
        )
        p=plot(trace, layout)
        display(p)
    end

    return (s, t)
end
foname="./StretchedDataSmallNoise.jld2"
# make Ricker wavelet

dt     = 0.05; # s
#f      = 0.05; # Hz
f = 1.0
w,tw = ricker(f=f, n=41, dt=dt, nargout=false);

## make random reflectivity sequence and convolve with Ricker
#npts = 601
npts = 6001; # length of time series

rng = MersenneTwister(20190603)
f    = randn(rng, Float64, npts);
u0   = conv(f,w)[1:npts]; # make waveform time series

## make time varying stretched + Noise
dvV = - 0.0001 #dv over V [%] = dt over T, therfore dt = dv/V * T

tvec  = collect(( 0:npts - 1 ) .* dt);
st   = tvec .* dvV
# add random noise (5 % of maximum st)
rng = MersenneTwister(20190603);
noise = randn(rng, Float64, npts);
noise = 0.05*maximum(abs.(st)).*(noise./maximum(abs.(st)))

#------------------#
#if adding noise on time shift
st = st .+ noise;
#------------------#

tvec2 = tvec .+ st;

spl = Spline1D(tvec, u0)
u1 = spl(tvec2);

st = st;

@save foname dt u0 u1 st

nargout = true

if nargout

    #time shift curve

    function lineplot1()
        tr = scatter(;x=tvec, y=st, mode="lines", name="Timeshift")
        layout = Layout(
            xaxis=attr(title="Time"),
            yaxis=attr(title="Timeshift"),
        )
        plot(tr, layout)
    end

    function lineplot2()
        tr1 = scatter(;x=tvec, y=u0, mode="lines+markers", name="u0")
        tr2 = scatter(;x=tvec, y=u1, mode="lines+markers", name="u1")
        layout = Layout(
            xaxis=attr(title="Time"),
            yaxis=attr(title="Signal"),
        )
        plot([tr1, tr2], layout)
    end

    p = [lineplot1(), lineplot2()]

end
