include("makeStretchData.jl")

using PlotlyJS, SeisIO, SeisNoise, JLD2

foname="verificationData.jld2"

finame_xcorr="../../xcorr_BP/testData/BPnetwork_Jan03_xcorrs.2003.1.T00:00:00.jld2"
dataset_xcorr="2003.1.T00:00:00/BP.CCRB..BP1.BP.EADB..BP1"
xcf=jldopen(finame_xcorr)
real_xcorr=xcf[dataset_xcorr]
save=false

if size(real_xcorr.corr)[2]>1 stack!(real_xcorr, allstack=true) end

dvVlist = collect(-0.03:0.0001:0.03)
noiselist = collect(0.0:0.001:0.1)

if save
    f=jldopen(foname, "a+")
    f["info/dvVlist"] = dvVlist
    f["info/noiselist"] = noiselist
end

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

        signal1_ds, t_ds = generateSignal("dampedSinusoid", params=dampedSinParams)
        signal2_ds, st_ds = stretchData(signal1_ds, dampedSinParams["dt"], dvV, n=noiselvl*10)

        signal1_ds = addNoise(signal1_ds, noiselvl)
        signal2_ds = addNoise(signal2_ds, noiselvl)

        # Example of ricker wavelet generation and convolution with random reflectivity
        # series, stretching, and noise addition.
        rickerParams = Dict( "f"     => 0.25,
                             "dt"    => 0.05,
                             "npr"   => 4001,
                             "npts"  => 4001)

        signal1_rc, t_rc = generateSignal("ricker", sparse=100, params=rickerParams)
        signal2_rc, st_rc = stretchData(signal1_rc, rickerParams["dt"], dvV, n=noiselvl*10)

        signal1_rc = addNoise(signal1_rc, noiselvl)
        signal2_rc = addNoise(signal2_rc, noiselvl)

        # Example of stretching real cross-correlations
        stretch_xcorr, st_xcorr = stretchData(real_xcorr.corr[:,1], 1/real_xcorr.fs, dvV, starttime=-100.0, stloc=0.0, n=noiselvl*10)

        if save
            f["dampedSinusoid/$dvV.$noiselvl"] = [signal1_ds, signal2_ds]
            f["rickerConv/$dvV.$noiselvl"] = [signal1_rc, signal2_rc]
            f["realData/$dvV.$noiselvl"] = [real_xcorr.corr[:,1], stretch_xcorr]
        end
    end
end
if save close(f) end
close(xcf)
