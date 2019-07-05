include("makeStretchData.jl")

using PlotlyJS, SeisIO, SeisNoise, JLD2

foname="verificationData.jld2"
f=jldopen(foname, "a+")

finame_xcorr="../../xcorr_BP/testData/BPnetwork_Jan03_xcorrs.2003.1.T00:00:00.jld2"
dataset_xcorr="2003.1.T00:00:00/BP.CCRB..BP1.BP.EADB..BP1"
xcf=jldopen(finame_xcorr)
real_xcorr=xcf[dataset_xcorr]

if size(real_xcorr.corr)[2]>1 stack!(real_xcorr, allstack=true) end

dvVlist = collect(-0.03:0.0001:0.03)
noiselist = collect(0.0:0.001:0.1)
f["info/dvVlist"] = dvVlist
f["info/noiselist"] = noiselist
for dvV in dvVlist
    for noiselvl in noiselist
        # Example of damped sinusoid generation, stretching, and noise addition
        dampedSinParams = Dict( "A"    => 1.0,
                                "ω"    => 0.25,
                                "ϕ"    => 0.0,
                                "λ"    => 0.02,
                                "dt"   => 0.05,
                                "η"    => 100.0,
                                "t0"   => 0.0,
                                "npts" => 4001)

        signal1, t = generateSignal("dampedSinusoid", params=dampedSinParams)
        signal2, st = stretchData(signal1, dampedSinParams["dt"], dvV, n=noiselvl*10)

        signal1 = addNoise(signal1, noiselvl)
        signal2 = addNoise(signal2, noiselvl)

        f["dampedSinusoid/$dvV.$noiselvl"] = [signal1, signal2]

        # Example of ricker wavelet generation and convolution with random reflectivity
        # series, stretching, and noise addition.
        rickerParams = Dict( "f"     => 0.25,
                             "dt"    => 0.05,
                             "npr"   => 4001,
                             "npts"  => 4001)

        signal1, t = generateSignal("ricker", sparse=100, params=rickerParams)
        signal2, st = stretchData(signal1, rickerParams["dt"], dvV, n=noiselvl*10)

        signal1 = addNoise(signal1, noiselvl)
        signal2 = addNoise(signal2, noiselvl)

        f["rickerConv/$dvV.$noiselvl"] = [signal1, signal2]

        # Example of stretching real cross-correlations
        stretch_xcorr, st = stretchData(real_xcorr.corr[:,1], 1/real_xcorr.fs, dvV, starttime=-100.0, stloc=0.0, n=noiselvl*10)

        f["realData/$dvV.$noiselvl"] = [real_xcorr.corr[:,1], stretch_xcorr]
    end
end
close(f)
close(xcf)
