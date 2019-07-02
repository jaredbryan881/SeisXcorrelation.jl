using SeisIO, SeisNoise, JLD2, Distributed, Dates

include("../../src/SeisXcorrelation.jl")
using SeisXcorrelation
include("../../src/pairing.jl")

# input parameters
InputDict = Dict( "finame"     => "./outputData/BPnetworkxcorr_Jan03.jld2",
                  "basefoname" => "testData",
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
data = jldopen(InputDict["finame"])
# read station and time stamp lists
stlist = data["info/stationlist"][:]
tstamplist = data["info/timestamplist"][1:2]

station_pairs = data["info/corrstationlist"]["xcorr"]
corrnames = station_pairs[1, :] .* "." .* station_pairs[2, :]
c3_pairs = generate_pairs(corrnames)
println(size(station_pairs))
println(station_pairs[:, 1])
exit()
# create output file and save station and pairing information in JLD2
jldopen("$(InputDict["basefoname"]).jld2", "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = station_pairs;
    file["info/tserrors"]        = []
end


for i=1:length(tstamplist)
    st = time()
    InputDict["foname"] = "$(InputDict["basefoname"])$i.jld2"
    # In this case, allowable_pairs is all of the precomputed xcorrs, which is all possible station pairs.
    # In this case, allowable_vs is all stations in the station list.
    errors = seisxcorrelation_highorder(tstamplist[i], station_pairs, corrnames, stlist, data, InputDict)

    jldopen("$(InputDict["basefoname"]).jld2", "a+") do file
        append!(file["info/tserrors"], errors)
    end
    et = time()
    println("$(tstamplist[i]) completed in $(et-st) seconds.")
end
close(data)
println("Successfully completed cross-correlation and saved to $(InputDict["foname"])")
#pmap(x -> seisxcorrelation(x, finame, foname, corrtype, corrorder, maxtimelag, freqmin, freqmax, fs, cc_len, cc_step), [tstamplist[1]])
