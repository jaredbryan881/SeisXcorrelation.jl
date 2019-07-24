using SeisIO, JLD2, Wavelets

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
    siglvl  = InputDict["siglvl"]
    twindow = [tmin, tmax]

    for dvV in dvVlist#[1:5:end]
        println("dv/v: $dvV")
        for noiselvl in noiselist
            signals = valData["$type/$dvV.$noiselvl"]
            ref = signals["ref"]
            cur = signals["cur"]

            npts = length(cur)
            tvec = (collect(0:npts-1) .* dt) .+ mintime

            (freqbands, dvv, err) = wxs_freqbands(cur, ref, tvec, twindow, fbands, dj, s0, J, sig=sig, wvn=wvn, unwrap=unwrap, siglvl=siglvl)
        end
    end
end

function test_WTS(finame::String, foname::String, InputDict::Dict, type::String, plot::Bool=true)
    println("test")
end

function test_WTDTW(finame::String, foname::String, InputDict::Dict, type::String; plot::Bool=true)
    println("test")
end

finame = "/Users/jared/SCECintern2019/data/validation/verificationData.jld2"
foname = "verificationDataError_WXS.jld2"
type = "rickerConv"

# setup parameters
wv = "ricker"
dvov = 0.005

InputDict_WXS = Dict( "fmin"    => 0.0,
                      "fmax"    => 10.0,
                      "tmin"    => 50.0,
                      "tmax"    => 100.0,
                      "dt"      => 0.05,
                      "dj"      => 1/12,
                      "s0"      => -1,
                      "J"       => -1,
                      "wvn"     => "morlet",
                      "sig"     => false,
                      "unwrap"  => true,
                      "mintime" => 0,
                      "siglvl"  => 0.95,
                      "fbands"  => [[0.1,0.2], [0.2, 0.4], [0.4, 0.8], [0.8, 1.6], [1.6, 3.2]])

for type in ["rickerConv", "dampedSin", "realData"]
    test_WXS(finame, foname, InputDict_WXS, type)
end
