using SeisIO, SeisNoise, JLD2, Distributed, Dates

@everywhere include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")

# input parameters
InputDict = Dict( "basefiname" => "./testData/BPnetwork_Jan03_xcorrs",
                  "basefoname" => "./testOutputc3/BPnetwork_Jan03_c3xcorrs",
                  "maxReadNum" => 6,
                  "start_lag"  => 0.0,
                  "window_len" => 100.0,
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 20,
                  "cc_step"    => 10,
                  "corrtype"   => ["xcorr"],
                  "corrorder"  => 3,
                  "maxtimelag" => 100.0,
                  "allstack"   => true)

# read data from JLD2
data = jldopen(InputDict["basefiname"]*".jld2")
# read station and time stamp lists
stlist = data["info/stationlist"]
tstamplist = data["info/timestamplist"]
station_pairs = data["info/corrstationlist"]["xcorr"]
corrnames = station_pairs[1, :] .* "." .* station_pairs[2, :]
c3_pairs = generate_pairs(corrnames)

# create output file and save station and pairing information in JLD2
jldopen("$(InputDict["basefoname"]).jld2", "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = station_pairs;
    file["info/tserrors"]        = []
end

# get indices of timesteps to load in each iteration
mapWindows = window_read(length(tstamplist), InputDict["maxReadNum"])
for win in mapWindows
    st=time()
    # load timesteps into array of dicts [Dict(ts/stn=>data for stn in stlist) for ts in tslist]
    fileHandles = [jldopen(InputDict["basefiname"]*".$i.jld2") for i in tstamplist[win]]
    dsets = read_JLD22Dict(fileHandles, tstamplist[win])

    # map seisxcorrelation over timesteps
    errors = pmap((x,y)->seisxcorrelation_highorder(x, y, station_pairs, corrnames, stlist, InputDict), dsets, tstamplist[win])

    for f in fileHandles close(f) end

    #TODO: write tserrors to basefile
    et=time()
    println("Section took $(et-st) seconds")
end

close(data)
println("Successfully completed cross-correlation and saved to $(InputDict["basefoname"])")
rmprocs(workers())
