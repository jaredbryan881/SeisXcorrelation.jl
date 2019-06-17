using SeisIO, Noise, JLD2, PlotlyJS

include("../../src/utils.jl")
include("../../src/Plotting/plot_data.jl")

# input parameters
InputDict = Dict( "corrname"   => "outputData/BPnetworkxcorr_neq.jld2",
                  "refname"    => "refXCorrs/reference_xcorr_neq.jld2",
                  "method"     => "MWCS",
                  "freqmin"    => 0.1,
                  "freqmax"    => 9.9,
                  "fs"         => 20.0,
                  "maxtimelag" => 100.0,
                  "mintimelag" => -100.0,
                  "dtt_lag"    => "static",
                  "dtt_v"      => 1.0,
                  "dtt_minlag" => -100.0,
                  "dtt_width"  => 10.0,
                  "window_step"=> 2.0,
                  "dtt_sides"  => "Both")

# read daily xcorrs and reference xcorr
xcorrs = jldopen(InputDict["corrname"])
reference = jldopen(InputDict["refname"])

tstamplist = xcorrs["info/timestamplist"][1:end-2]
stationpairs = reference["info/reference_xcorrnames"]

#iterate over time steps
for tstamp in tstamplist
    for stnpair in stationpairs[37:end]
        ref = reference[stnpair]
        cur = xcorrs["$tstamp/$stnpair"]
        println("$tstamp/$stnpair")
        stack!(cur, allstack=true)

        dist = xcorrs["$tstamp/$stnpair"].misc["dist"]
        println(xcorrs["$tstamp/$stnpair"].loc)
        data = jldopen("/Users/jared/SCECintern2019/RemoveEarthquakes/dataset/BPnetwork_RemovedEQ.jld2")
        println(data["$tstamp/BP.CCRB..BP1"].loc)
        println(data["$tstamp/BP.EADB..BP1"].loc)
        exit()

        if InputDict["method"] == "MWCS"
            time_axis, dt, err, coh = mwcs(ref/maximum(ref),
                                           cur.corr/maximum(cur.corr),
                                           InputDict["freqmin"],
                                           InputDict["freqmax"],
                                           InputDict["fs"],
                                           InputDict["mintimelag"],
                                           InputDict["dtt_width"],
                                           InputDict["window_step"],
                                           0)
            m, em, a, ea, m0, em0 = mwcs_dvv(time_axis,
                                             dt,
                                             err,
                                             coh,
                                             InputDict["dtt_lag"],
                                             dist/1000,
                                             InputDict["dtt_v"],
                                             InputDict["dtt_minlag"],
                                             InputDict["dtt_width"],
                                             InputDict["dtt_sides"])

        elseif InputDict["method"] == "stretching"
            println("Using Stretching")
            dv, cc, cdp, eps, err, C = stretching(ref, cur, time_axis, window, InputDict["freqmin"], InputDict["freqmax"], -0.03, 0.03, 100)
        end

        # save dv/v to a GeoDataFrame (?)
    end
end

close(xcorrs)
close(ref)
