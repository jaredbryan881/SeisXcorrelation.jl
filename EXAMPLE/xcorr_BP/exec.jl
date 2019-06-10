include("../../src/SeisXcorrelation.jl")
include("../../src/pairing.jl")

using SeisIO, Noise, JLD2, Distributed

# input parameters
finame = "../../../SeisDownload.jl/EXAMPLE/Download_BP/dataset/BPnetwork.jld2"
foname = "dataOutput.jld2"
corrorder = 1
corrtype = ["xcorr", "xchancorr", "acorr"]
maxtimelag = 100.0
freqmin = 0.1
freqmax = 10.0
fs = 20.0
cc_len = 3600
cc_step = 1800

# read data from JLD2
data = jldopen(finame)

# read station and time stamp lists
stlist = data["info/stationlist"]
tstamplist = data["info/timestamplist"]

# generate station pairs
if corrorder == 1
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
jldopen(foname, "w") do file
    file["info/timestamplist"]   = tstamplist;
    file["info/stationlist"]     = stlist;
    file["info/corrstationlist"] = sorted_pairs;
end

pmap(x -> seisxcorrelation(x, finame, foname, corrtype, corrorder, maxtimelag, freqmin, freqmax, fs, cc_len, cc_step), tstamplist)
