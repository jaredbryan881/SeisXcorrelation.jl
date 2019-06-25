using SeisIO, Noise, JLD2, Distributed, Dates

include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")

# input parameters
InputDict = Dict( "finame"     => "./dataset/BPnetwork_Jan03.jld2",
                  "basefoname" => "testData",
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "cc_len"     => 3600,
                  "cc_step"    => 1800,
                  "corrtype"   => ["xcorr", "acorr", "xchancorr"],
                  "corrorder"  => 1,
                  "maxtimelag" => 100.0,
                  "allstack"   => true)

# read data from JLD2
data = jldopen(InputDict["finame"])
# read station and time stamp lists
stlist = data["info/stationlist"][:]
tstamplist = data["info/DLtimestamplist"][1:2]

station_pairs = generate_pairs(stlist)
# sort station pairs into autocorr, xcorr, and xchancorr
sorted_pairs = sort_pairs(station_pairs)

# create output file and save station and pairing information in JLD2
jldopen("$(InputDict["basefoname"])*.jld2", "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = sorted_pairs;
    file["info/tserrors"]        = []
end

# TODO make sure tserrors are actually written to file
for i=1:length(tstamplist)
    st = time()
    InputDict["foname"] = "$(InputDict["basefoname"])$i.jld2"
    errors = seisxcorrelation(tstamplist[i], stlist, data, InputDict)

    jldopen("$(InputDict["basefoname"]).jld2", "a+") do file
        append!(file["info/tserrors"], errors)
    end

    et = time()
    println("$(tstamplist[i]) completed in $(et-st) seconds.")
end
close(data)
println("Successfully completed cross-correlation and saved to $(InputDict["foname"])")
#pmap(x -> seisxcorrelation(x, finame, foname, corrtype, corrorder, maxtimelag, freqmin, freqmax, fs, cc_len, cc_step), [tstamplist[1]])
