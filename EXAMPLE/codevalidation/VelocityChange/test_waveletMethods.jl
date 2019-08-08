using SeisIO, JLD2, Wavelets, Statistics, DSP, FFTW, Plots, GLM, DataFrames

include("wavelet_methods.jl")

function test_WXS(finame::String, foname::String, InputDict::Dict, type::String; plot::Bool=true)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    fmin    = InputDict["fmin"]
    fmax    = InputDict["fmax"]
    tmin    = InputDict["tmin"]
    tmax    = InputDict["tmax"]
    dt      = InputDict["dt"]
    dj      = InputDict["dj"]
    s0      = InputDict["s0"]
    J       = InputDict["J"]
    wvn     = InputDict["wvn"]
    sig     = InputDict["sig"]
    fbands  = InputDict["fbands"]
    unwrap  = InputDict["unwrap"]
    mintime = InputDict["mintime"]
    winflag = InputDict["winflag"]
    nwindow = InputDict["nwindow"]
    siglvl  = InputDict["siglvl"]
    twindow = [tmin, tmax]
    fwindow = [fmin, fmax]

    # output file
    f_out = jldopen(foname, "w")

    for dvV in dvVlist
        println("dv/v: $dvV")
        for noiselvl in noiselist
            signals = valData["$type/$dvV.$noiselvl"]
            ref = signals["ref"]
            cur = signals["cur"]

            npts = length(cur)
            tvec = (collect(0:npts-1) .* dt) .+ mintime

            # Wavelet cross-spectrum on discrete frequency bands
            (freqbands, dvv, err) = wxs_freqbands(cur, ref, tvec, twindow, fbands, dj, s0, J, wvn=wvn, unwrapflag=unwrap, sig=sig, siglvl=siglvl, windowflag=winflag, nwindow=nwindow)
            f_out["bands/$dvV/$noiselvl"] = [freqbands, dvv, err]
            p=PlotlyJS.plot(dvv, mode="lines+markers")
            display(p)
            readline()

            # Wavelet cross-spectrum for each frequency
            (freqs, dvv, err) = wxs_allfreq(cur, ref, tvec, twindow, fwindow, dj, s0, J, wvn=wvn, unwrapflag=unwrap, sig=sig, siglvl=siglvl, windowflag=winflag, nwindow=nwindow)
            f_out["full/$dvV/$noiselvl"] = [freqbands, dvv, err]
            p=PlotlyJS.plot(freqs, dvv, mode="lines+markers")
            display(p)
            readline()
        end
    end
    close(f_out)
end

function test_WTS(finame::String, foname::String, InputDict::Dict, type::String, plot::Bool=true)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    fmin    = InputDict["fmin"]
    fmax    = InputDict["fmax"]
    tmin    = InputDict["tmin"]
    tmax    = InputDict["tmax"]
    dt      = InputDict["dt"]
    dj      = InputDict["dj"]
    s0      = InputDict["s0"]
    J       = InputDict["J"]
    wvn     = InputDict["wvn"]
    ndv     = InputDict["ndv"]
    norm    = InputDict["norm"]
    dvmax   = InputDict["dvmax"]
    fbands  = InputDict["fbands"]
    mintime = InputDict["mintime"]
    winflag = InputDict["winflag"]
    nwindow = InputDict["nwindow"]

    twindow = [tmin, tmax]
    fwindow = [fmin, fmax]

    # output file
    f_out = jldopen(foname, "w")

    for dvV in dvVlist
        println("dv/v: $dvV")
        for noiselvl in noiselist
            signals = valData["$type/$dvV.$noiselvl"]
            ref = signals["ref"]
            cur = signals["cur"]

            npts = length(cur)
            tvec = (collect(0:npts-1) .* dt) .+ mintime

            # Wavelet trace stretching on discrete frequency bands
            (freqbands, dvv, cc, cdp, err) = wts_freqbands(cur, ref, tvec, twindow, fbands, dj, s0, J, dvmax=dvmax, ndv=ndv, wvn=wvn, normalize=norm, windowflag=winflag, nwindow=nwindow)
            f_out["bands/$dvV/$noiselvl"] = [freqbands, dvv, err, cc, cdp]

            # Wavelet trace stretching for each frequency
            (freqs, dvv, cc, cdp, err) = wts_allfreqs(cur, ref, tvec, twindow, fbands, dj, s0, J, dvmax=dvmax, ndv=ndv, wvn=wvn, normalize=norm, windowflag=winflag, nwindow=nwindow)
            f_out["full/$dvV/$noiselvl"] = [freqs, dvv, err, cc, cdp]
        end
    end
    close(f_out)
end

function test_WTDTW(finame::String, foname::String, InputDict::Dict, type::String; plot::Bool=true)
    valData = jldopen(finame)
    dvVlist = valData["info/dvVlist"]
    noiselist = valData["info/noiselist"]

    fmin    = InputDict["fmin"]
    fmax    = InputDict["fmax"]
    tmin    = InputDict["tmin"]
    tmax    = InputDict["tmax"]
    dt      = InputDict["dt"]
    dj      = InputDict["dj"]
    s0      = InputDict["s0"]
    J       = InputDict["J"]
    b       = InputDict["b"]
    wvn     = InputDict["wvn"]
    ndv     = InputDict["ndv"]
    dir     = InputDict["dir"]
    norm    = InputDict["norm"]
    dvmax   = InputDict["dvmax"]
    fbands  = InputDict["fbands"]
    maxLag  = InputDict["maxLag"]
    winflag = InputDict["winflag"]
    mintime = InputDict["mintime"]
    nwindow = InputDict["nwindow"]

    twindow = [tmin, tmax]
    fwindow = [fmin, fmax]

    # output file
    f_out = jldopen(foname, "w")

    for dvV in dvVlist
        println("dv/v: $dvV")
        for noiselvl in noiselist
            signals = valData["$type/$dvV.$noiselvl"]
            ref = signals["ref"]
            cur = signals["cur"]

            npts = length(cur)
            tvec = (collect(0:npts-1) .* dt) .+ mintime

            # Wavelet dynamic time warping on discrete frequency bands
            (freqbands, dvv, err) = wtdtw_freqbands(cur, ref, tvec, twindow, fbands, maxLag, b, dir, dj, s0, J, wvn=wvn, normalize=norm, windowflag=winflag, nwindow=nwindow)
            f_out["bands/$dvV/$noiselvl"] = [freqbands, dvv, err]

            # Wavelet dynamic time warping for each frequency
            (freqs, dvv, err) = wtdtw_allfreqs(cur, ref, tvec, twindow, fwindow, maxLag, b, dir, dj, s0, J, wvn=wvn, normalize=norm, windowflag=winflag, nwindow=nwindow)
            f_out["full/$dvV/$noiselvl"] = [freqs, dvv, err]
        end
    end
    close(f_out)
end

finame = "/Users/jared/SCECintern2019/data/validation/verificationData.jld2"
foname_WXS = "verificationData_WXS.jld2"
foname_WTS = "verificationData_WTS.jld2"
foname_DTW = "verificationData_DTW.jld2"

InputDict_WXS = Dict( "fmin"    => 0.0,
                      "fmax"    => 10.0,
                      "tmin"    => 0.0,
                      "tmax"    => 200.0,
                      "dt"      => 0.05,
                      "dj"      => 1/12,
                      "s0"      => -1,
                      "J"       => -1,
                      "wvn"     => "morlet",
                      "sig"     => false,
                      "unwrap"  => false,
                      "mintime" => 0,
                      "siglvl"  => 0.95,
                      "fbands"  => [0.1 0.2; 0.2 0.4; 0.4 0.8; 0.8 1.6; 1.6 3.2],
                      "winflag" => true,
                      "nwindow" => 1.5)

InputDict_WTS = Dict( "fmin"    => 0.0,
                      "fmax"    => 10.0,
                      "tmin"    => 0.0,
                      "tmax"    => 200.0,
                      "dt"      => 0.05,
                      "dj"      => 1/12,
                      "s0"      => -1,
                      "J"       => -1,
                      "wvn"     => "morlet",
                      "ndv"     => 100,
                      "dvmax"   => 0.1,
                      "mintime" => 0,
                      "norm"    => true,
                      "fbands"  => [0.1 0.2; 0.2 0.4; 0.4 0.8; 0.8 1.6; 1.6 3.2],
                      "winflag" => true,
                      "nwindow" => 1.5)

InputDict_DTW = Dict( "fmin"    => 0.0,
                      "fmax"    => 10.0,
                      "tmin"    => 0.0,
                      "tmax"    => 200.0,
                      "dt"      => 0.05,
                      "dj"      => 1/12,
                      "s0"      => -1,
                      "b"       => 1,
                      "J"       => -1,
                      "wvn"     => "morlet",
                      "dir"     => 0,
                      "ndv"     => 100,
                      "dvmax"   => 0.1,
                      "mintime" => 0,
                      "norm"    => true,
                      "maxLag"  => 150,
                      "fbands"  => [0.1 0.2; 0.2 0.4; 0.4 0.8; 0.8 1.6; 1.6 3.2],
                      "winflag" => true,
                      "nwindow" => 1.5)


for type in ["spectSynth", "rickerConv", "dampedSin", "realData"]
    test_WXS(finame, foname_WXS, InputDict_WXS, type)
    test_WTS(finame, foname_WTS, InputDict_WTS, type)
    test_WTDTW(finame, foname_DTW, InputDict_DTW, type)
end
