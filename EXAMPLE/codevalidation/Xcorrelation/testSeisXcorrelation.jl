using SeisIO, JLD2, SeisNoise, FFTW, Dates

include("../../../src/SeisXcorrelation.jl")
include("../../../src/pairing.jl")

# This is a simple test script for basic Noise.jl functionality.
# This is the core process used in SeisXcorrelation.jl
# Read data from JLD2 -> Package as SeisData -> FFT -> XCorr -> Save as JLD2

InputDict = Dict( "finame"     => "./sinStretchTest.jld2",
                  "basefoname" => "testData",
                  "timeunit"   => 60,
                  "freqmin"    => 0.01,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 58,
                  "cc_step"    => 1,
                  "corrtype"   => ["xcorr", "acorr", "xchancorr"],
                  "corrorder"  => 1,
                  "maxtimelag" => 100.0,
                  "allstack"   => true)

data=jldopen(InputDict["finame"])
stlist = data["info/stationlist"]
tstamplist = data["info/DLtimestamplist"]
station_pairs = generate_pairs(stlist)
sorted_pairs = sort_pairs(station_pairs)

# create output file and save station and pairing information in JLD2
jldopen((InputDict["basefoname"])*".jld2", "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = sorted_pairs;
    file["info/tserrors"]        = []
end

dsets = read_JLD22Dict(data, tstamplist)
# compute xcorr
errors = seisxcorrelation(dsets[1], tstamplist[1], InputDict)
close(data)
