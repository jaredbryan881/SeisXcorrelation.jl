using SeisIO, FileIO, JLD2, SeisNoise, FFTW, PlotlyJS, Dates

include("../../../src/SeisXcorrelation.jl")
include("../../../src/pairing.jl")

# This is a simple test script for basic Noise.jl functionality.
# This is the core process used in SeisXcorrelation.jl
# Read data from JLD2 -> Package as SeisData -> FFT -> XCorr -> Save as JLD2

# read time series data from JLD2
finame="sinStretchTest.jld2"

InputDict = Dict( "finame"     => "./sinStretchTest.jld2",
                  "foname"     => "testData.jld2",
                  "freqmin"    => 0.01,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 58,
                  "cc_step"    => 1,
                  "corrtype"   => ["xcorr", "acorr", "xchancorr"],
                  "corrorder"  => 1,
                  "maxtimelag" => 100.0,
                  "allstack"   => true)

data=jldopen(finame)
stlist = data["info/stationlist"]
tstamplist = data["info/DLtimestamplist"]
station_pairs = generate_pairs(stlist)
sorted_pairs = sort_pairs(station_pairs)

# create output file and save station and pairing information in JLD2
jldopen((InputDict["foname"]), "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = sorted_pairs;
    file["info/tserrors"]        = []
end

# compute xcorr
errors = seisxcorrelation(tstamplist[1], stlist, data, InputDict)
close(data)
