include("makeStretchData.jl")
include("../../../src/utils.jl")

using PlotlyJS, SeisIO, SeisNoise, JLD2, FFTW, ORCA

foname="/Users/jared/SCECintern2019/data/validation/verificationData.jld2"

types = ["dampedSinusoid", "rickerConv", "realData", "spectSynth"]

finame_xcorr="/Users/jared/SCECintern2019/data/reference/BPnetwork_200345_ref.jld2"
dataset_xcorr="BP.CCRB..BP1.BP.SMNB..BP1"
xcf=jldopen(finame_xcorr)
real_xcorr=xcf[dataset_xcorr]
save=true
plot=false

if size(real_xcorr.corr)[2]>1 stack!(real_xcorr, allstack=true) end

dvVlist = collect(-0.03:0.005:0.03)
noiselist = collect(0:0.01:0.1)

if save
    f=jldopen(foname, "a+")
    f["info/dvVlist"] = dvVlist
    f["info/noiselist"] = noiselist
end

dampedSinParams = Dict( "A"    => 1.0,
                        "ω"    => 0.75,
                        "ϕ"    => 0.0,
                        "λ"    => 0.025,
                        "dt"   => 0.05,
                        "η"    => 0.0,
                        "t0"   => 0.0,
                        "npts" => 4001 )

sincParams = Dict( "A"    => 1.0,
                   "ω"    => 0.1,
                   "ϕ"    => 0.0,
                   "dt"   => 0.05,
                   "t0"   => 0.0,
                   "npts" => 4001 )

rickerParams = Dict( "f"      => 0.25,
                     "dt"     => 0.05,
                     "npr"    => 4001,
                     "npts"   => 4001,
                     "m"      => 1,
                     "sparse" => 0,
                     "stretchSource" => 0.0 )

chirpParams = Dict( "c"     => collect(range(0.15, stop=15.0, length=100)),
                    "tp"    => collect(range(0.04, stop=4.0, length=100)),
                    "mintp" => 0.04,
                    "maxtp" => 4.0,
                    "dist"  => 1000.0,
                    "n"     => 500,
                    "dt"    => 0.002,
                    "t0"    => 0.0 )

spectParams = Dict( "A"    => 1.0,
                    "dt"   => 0.05,
                    "npts" => 4001,
                    "t0"   => 0.0 )

lags = -real_xcorr.maxlag:1/real_xcorr.fs:real_xcorr.maxlag

for dvV in dvVlist
    for noiselvl in noiselist
        if "dampedSinusoid" in types
            # Example of damped sinusoid generation, stretching, and noise addition
            signal1_ds, t_ds = generateSignal("dampedSinusoid", dampedSinParams)
            normalize!(signal1_ds)
            signal2_ds, st_ds = stretchData(signal1_ds, dampedSinParams["dt"], dvV, n=noiselvl)

            addNoise!(signal1_ds, noiselvl)
            addNoise!(signal2_ds, noiselvl, seed=664739)
        end

        if "rickerConv" in types
            # Example of ricker wavelet generation and convolution with random reflectivity
            # series, stretching, and noise addition
            signal1_rc, t_rc = generateSignal("ricker", rickerParams, sparse=rickerParams["sparse"], stretchSource=0.0)
            normalize!(signal1_rc)
            if rickerParams["stretchSource"] != 0.0
                inSignal, t_rc = generateSignal("ricker", rickerParams, sparse=rickerParams["sparse"], stretchSource=rickerParams["stretchSource"])
                normalize!(inSignal)
            else
                inSignal = signal1_rc
            end
            signal2_rc, st_rc = stretchData(inSignal, rickerParams["dt"], dvV, n=noiselvl)

            addNoise!(signal1_rc, noiselvl)
            addNoise!(signal2_rc, noiselvl, seed=664739)
        end

        if "spectSynth" in types
            # Example of synthetic generation in frequency domain and irfft to get into time domain
            signal1_ss, t_ss = generateSignal("spectSynth", spectParams)
            normalize!(signal1_ss)
            signal2_ss, st_ss = stretchData(signal1_ss, spectParams["dt"], dvV, n=noiselvl)

            addNoise!(signal1_ss, noiselvl)
            addNoise!(signal2_ss, noiselvl, seed=664739)
        end

        if "realData" in types
            # Example of stretching real cross-correlations
            xcorr = real_xcorr.corr[:,1]
            normalize!(xcorr)
            stretch_xcorr, st_xcorr = stretchData(xcorr, 1/real_xcorr.fs, dvV, starttime=-real_xcorr.maxlag, stloc=0.0, n=noiselvl)

            addNoise!(xcorr, noiselvl)
            addNoise!(stretch_xcorr, noiselvl, seed=664739)
        end

        if save
            # damped sinusoid
            f["dampedSinusoid/$dvV.$noiselvl/ref"] = signal1_ds
            f["dampedSinusoid/$dvV.$noiselvl/cur"] = signal2_ds
            # ricker conv
            f["rickerConv/$dvV.$noiselvl/ref"] = signal1_rc
            f["rickerConv/$dvV.$noiselvl/cur"] = signal2_rc
            # real data
            f["realData/$dvV.$noiselvl/ref"] = xcorr
            f["realData/$dvV.$noiselvl/cur"] = stretch_xcorr
            # spect synth
            f["spectSynth/$dvV.$noiselvl/ref"] = signal1_ss
            f["spectSynth/$dvV.$noiselvl/cur"] = signal2_ss
        end

        if plot
            plots = Array{PlotlyJS.Plot, 1}()
            if "dampedSinusoid" in types
                p1 = PlotlyJS.Plot([PlotlyJS.scatter(;x=t_ds, y=signal1_ds, name="Unstretched", line_color="black"),
                                    PlotlyJS.scatter(;x=t_ds, y=signal2_ds, name="Stretched $(dvV*(-100))%", line_color="red")])
            else
                p1=PlotlyJS.Plot([PlotlyJS.scatter(;y=ones(10))])
            end
            if "rickerConv" in types
                p2 = PlotlyJS.Plot([PlotlyJS.scatter(;x=t_rc, y=signal1_rc, name="Unstretched", line_color="black"),
                                    PlotlyJS.scatter(;x=t_rc, y=signal2_rc, name="Stretched $(dvV*(-100))%", line_color="red")])
            else
                p2=PlotlyJS.Plot([PlotlyJS.scatter(;y=ones(10))])
            end
            if "realData" in types
                p3 = PlotlyJS.Plot([PlotlyJS.scatter(;x=lags, y=xcorr, name="Unstretched", line_color="black"),
                                    PlotlyJS.scatter(;x=lags, y=stretch_xcorr, name="Stretched $(dvV*(-100))%", line_color="red")])
            else
                p3=PlotlyJS.Plot([PlotlyJS.scatter(;y=ones(10))])
            end
            # temporarily, you cannot plot only certain pairs of data types. You must plot all.
            # these conditionals will be more meaningful later
            plots=[p3]

            p=PlotlyJS.plot(plots)
            display(p)
            readline()
        end
    end
end
if save close(f) end
close(xcf)
