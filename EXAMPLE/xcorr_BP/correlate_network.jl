using SeisIO, SeisNoise, JLD2, Distributed, Dates

@everywhere include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")
include("../../src/utils.jl")

# input parameters
InputDict = Dict( "finame"     => "/Users/jared/SCECintern2019/data/seismic/BPnetwork_Jan03_staltakurt_1.5_3.0.jld2",
                  "timeunit"   => 86400,     # unit of time xcorrs are saved in. DL_time_unit if SeisDownload was used
                  "basefoname" => "/Users/jared/SCECintern2019/data/xcorrs/BPnetwork_Jan03_xcorrs_1.5_3.0",
                  "maxReadNum" => 9,         # number of processes to read from JLD2 to dict and map over
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.99,
                  "fs"         => 20.0,
                  "cc_len"     => 3600,
                  "cc_step"    => 1800,
                  "corrtype"   => ["acorr", "xchancorr", "xcorr"],
                  "maxtimelag" => 100.0,
                  "to_whiten"  => true,      # true, false
                  "time_norm"  => "one-bit", # false, phase, or one-bit
                  "allstack"   => false)     # for now, only stacking over DL_time_unit is supported pre-save

# read data from JLD2
data = jldopen(InputDict["finame"])
# read station and time stamp lists
stlist = data["info/stationlist"]
tstamplist = data["info/DLtimestamplist"]

station_pairs = generate_pairs(stlist) # array
# sort station pairs into autocorr, xcorr, and xchancorr
sorted_pairs = sort_pairs(station_pairs) # dict

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
