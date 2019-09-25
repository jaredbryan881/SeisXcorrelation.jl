using SeisIO, SeisNoise, JLD2, Distributed, Dates

@everywhere include("../../src/SeisXcorrelation.jl")
#@everywhere using SeisXcorrelation
include("../../src/pairing.jl")
include("../../src/utils.jl")
include("../../src/reference.jl")
include("../../src/stacking.jl")
# input parameters
InputDict = Dict( "finame"     => "/Volumes/Elements/dataSave/BPnetwork/BPnetwork_2003/BPnetwork_2003_neq.jld2",
                  "timeunit"   => 86400,     # unit of time xcorrs are saved in. DL_time_unit if SeisDownload was used
                  "basefoname" => "./data/test",
                  "maxReadNum" => 3,         # number of processes to read from JLD2->dict and map over
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.99,
                  "fs"         => 20.0,
                  "cc_len"     => 3600,
                  "cc_step"    => 600,
                  "corrtype"   => ["acorr", "xchancorr", "xcorr"],
                  "corrmethod" => "cross-correlation", # "coherence", "cross-correlation", or "deconv"
                  "half_win"   => 20,
                  "maxtimelag" => 100.0,
                  "to_whiten"  => false, # true, false
                  "time_norm"  => false, # false, "phase", or "one-bit"
                  "allstack"   => false) # for now, only linear stacking over DL_time_unit is supported pre-save

# read data from JLD2
data = jldopen(InputDict["finame"])
# read station and time stamp lists
stlist = data["info/stationlist"]
stlist = sort(stlist)
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
