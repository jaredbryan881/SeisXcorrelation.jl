using SeisIO, SeisNoise, JLD2, Distributed, Dates

@everywhere include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")

# input parameters
InputDict = Dict( "basefiname" => "/Users/jared/SCECintern2019/data/xcorrs/BPnetwork_Jan03_nowhiten_nonorm/BPnetwork_Jan03_xcorrs_1.5_3.0",
                  "basefoname" => "/Users/jared/SCECintern2019/data/xcorrs/BPnetwork_Jan03_c3/BPnetwork_Jan03_c3xcorrs",
                  "maxReadNum" => 6,
                  "start_lag"  => 0.0,
                  "window_len" => 100.0,
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 20,
                  "cc_step"    => 10,
                  "corrtype"   => ["xcorr"],
                  "maxtimelag" => 100.0,
                  "allstack"   => true)

# read data from JLD2
data = jldopen(InputDict["basefiname"]*".jld2")
# read station and time stamp lists
stlist = data["info/stationlist"]
tstamplist = data["info/timestamplist"]
station_pairs = data["info/corrstationlist"]["xcorr"]
close(data)
corrnames = station_pairs[1, :] .* "." .* station_pairs[2, :]

# create output file and save station and pairing information in JLD2
jldopen("$(InputDict["basefoname"]).jld2", "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = station_pairs;
end

# get indices of timesteps to load in each iteration
mapWindows = window_read(length(tstamplist), InputDict["maxReadNum"])
for win in mapWindows
    st=time()
    # load timesteps into array of dicts [Dict(ts/stn=>data for stn in stlist) for ts in tslist]
    fileHandles = [jldopen(InputDict["basefiname"]*".$i.jld2") for i in tstamplist[win]]
    dsets = read_JLD22Dict(fileHandles, tstamplist[win])

    # map seisxcorrelation over timesteps
    pmap((x,y)->seisxcorrelation_highorder(x, y, station_pairs, corrnames, stlist, InputDict), dsets, tstamplist[win])

    for f in fileHandles close(f) end

    et=time()
    println("Section took $(et-st) seconds")
end

println("Successfully completed cross-correlation and saved to $(InputDict["basefoname"])")
rmprocs(workers())
