using SeisIO, FileIO, JLD2, Noise, FFTW, PlotlyJS, Dates

# This is a simple test script for basic Noise.jl functionality. 
# This is the core process used in SeisXcorrelation.jl
# Read data from JLD2 -> Package as SeisData -> FFT -> XCorr -> Save as JLD2 

# read time series data from JLD2
finame="sinStretchTest.jld2"
f=jldopen(finame)
u1 = f["u1"]
u2 = f["u2"]

# fft params
freqmin = 0.01
freqmax = 9.9
fs      = 20.0
cc_len  = 58
cc_step = 1
# xcorr params
maxlag  = 100.0
N       = 1200

# put time series data into SeisChannel objects
# first SeisChannel
Ch1 = SeisChannel()
Ch1.x = u1
Ch1.name = "u1"
t1=Array{Int64, 2}(undef, 2, 2)
t1[1,1]=1
t1[1,2]=0 # arbitrary choice
t1[2,1]=length(u1)
t1[2,2]=0
Ch1.t=t1
Ch1.fs=fs
# second SeisChannel
Ch2 = SeisChannel()
Ch2.x = u2
Ch2.name = "u2"
t2=Array{Int64, 2}(undef, 2, 2)
t2[1,1]=1
t2[1,2]=0 # arbitrary choice
t2[2,1]=length(u2)
t2[2,2]=0
Ch2.t=t2
Ch2.fs=fs

# put SeisChannel objects into SeisData object
SD1 = SeisData(Ch1)
SD2 = SeisData(Ch2)

fft1 = compute_fft(SD1, freqmin, freqmax, fs, cc_step, cc_len)
fft2 = compute_fft(SD2, freqmin, freqmax, fs, cc_step, cc_len)
fft1.name = "fft1"
fft2.name = "fft2"

xcorr = compute_cc(fft2, fft1, maxlag)
tRange = collect(-length(xcorr.corr)/2:length(xcorr.corr)/2 - 1)*(1/fs)

@save "./xcorrValidation.jld2"
