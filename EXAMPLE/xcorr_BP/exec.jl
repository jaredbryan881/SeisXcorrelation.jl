using SeisIO, Noise, JLD2, Distributed, Dates

include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")

# input parameters
#finame = "../../../SeisDownload.jl/EXAMPLE/Download_BP/dataset/BPnetwork.jld2"

InputDict = Dict( "finame"     => "/Users/jared/SCECintern2019/RemoveEarthquakes/inputdata/BPnetwork.jld2",
                  "foname"     => "./outputData/BPnetworkxcorr_weq.jld2",
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 3600,
                  "cc_step"    => 1800,
                  "corrtype"   => ["xcorr", "xchancorr", "acorr"],
                  "corrorder"  => 1,
                  "maxtimelag" => 100.0,
                  "allstack"   => false)

# read data from JLD2
data = jldopen(InputDict["finame"])

# read station and time stamp lists
stlist = data["info/stationlist"][:]
tstamplist = data["info/DLtimestamplist"][:]

# generate station pairs
if InputDict["corrorder"] == 1
    station_pairs = generate_pairs(stlist)
    # sort station pairs into autocorr, xcorr, and xchancorr
    sorted_pairs = sort_pairs(station_pairs)
else
    station_pairs = data["info/corrstationlist"]["xcorr"]
    # concat correlated station names
    station_pairs = station_pairs[1, :] .* '.' .* station_pairs[2, :]
    corrstationlist = generate_pairs(station_pairs)
    # sort correlated pairs into autocorr, xcorr, and xchancorr
    # c2 and c3 will iterate only over xcorr
    sorted_pairs = sort_pairs(corrstationlist)
end

# create output file and save station and pairing information in JLD2
jldopen(InputDict["foname"], "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = sorted_pairs;
    file["info/tserrors"]        = []
end

for i=1:length(tstamplist)
    errors = seisxcorrelation(tstamplist[i], InputDict)

    jldopen(InputDict["foname"], "r+") do file
        append!(file["info/tserrors"], errors)
    end
end
close(data)
#pmap(x -> seisxcorrelation(x, finame, foname, corrtype, corrorder, maxtimelag, freqmin, freqmax, fs, cc_len, cc_step), [tstamplist[1]])
