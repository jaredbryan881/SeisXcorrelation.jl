using SeisIO, SeisNoise, JLD2, Distributed, Dates

@everywhere include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")
include("../../src/utils.jl")

# input parameters
InputDict = Dict( "finame"     => "/Users/jared/SCECintern2019/data/seismic/BPnetwork_Jan03.jld2",
                  "timeunit"   => 86400,
                  "basefoname" => "/Users/jared/SCECintern2019/data/xcorrs/BPnetwork_Jan03_weq_xcorrs",
                  "maxReadNum" => 8,
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 3600,
                  "cc_step"    => 1800,
                  "corrtype"   => ["acorr", "xchancorr", "xcorr"],
                  "corrorder"  => 1,
                  "maxtimelag" => 100.0,
                  "allstack"   => true)

# read data from JLD2
data = jldopen(InputDict["finame"])
# read station and time stamp lists
stlist = data["info/stationlist"]
tstamplist = data["info/DLtimestamplist"]

station_pairs = generate_pairs(stlist)
# sort station pairs into autocorr, xcorr, and xchancorr
sorted_pairs = sort_pairs(station_pairs)

# create base output file with network, station, and pairing info
jldopen("$(InputDict["basefoname"]).jld2", "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = sorted_pairs;
end

st=time()
# get indices of timesteps to load in each iteration
mapWindows = window_read(length(tstamplist), InputDict["maxReadNum"])
for win in mapWindows
    # load timesteps into array of dicts [Dict(ts/stn=>data for stn in stlist) for ts in tslist]
    dsets = read_JLD22Dict(data, tstamplist[win])
    # map seisxcorrelation over timesteps
    pmap((x,y)->seisxcorrelation(x, y, InputDict), dsets, tstamplist[win])
end
et=time()
close(data)
println("Successfully completed and saved cross-correlations in $(et-st) seconds.")
rmprocs(workers())
