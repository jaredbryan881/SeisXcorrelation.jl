using SeisIO, Noise, JLD2, PlotlyJS

include("../../src/utils.jl")

# input parameters
InputDict = Dict( "corrname"           => "outputData/BPnetworkxcorr_neq.jld2",
                  "refname"            => "refXCorrs/reference_xcorr_neq.jld2",
                  "method"             => "MWCS",
                  "freqmin"            => 0.1,
                  "freqmax"            => 9.9,
                  "fs"                 => 20.0,
                  "maxtimelag"         => 100.0,
                  "mintimelag"         => -100.0,
                  "dtt_lag"            => "static",
                  "dtt_v"              => 1.0,
                  "dtt_minlag"         => -100.0,
                  "dtt_width"          => 10.0,
                  "window_length"      => 10.0,
                  "window_step"        => 5.0,
                  "dtt_sides"          => "Both",
                  "smoothing_half_win" => 5)

# read daily xcorrs and reference xcorr
xcorrs = jldopen(InputDict["corrname"])
reference = jldopen(InputDict["refname"])

tstamplist = xcorrs["info/timestamplist"][:]
stationpairs = reference["info/reference_xcorrnames"]

#iterate over time steps
for tstamp in tstamplist[1:3]
    for stnpair in stationpairs[37:end]
        ref = reference[stnpair]
        cur = xcorrs["$tstamp/$stnpair"]
        println("$tstamp/$stnpair")
        stack!(cur, allstack=true)

        dist = xcorrs["$tstamp/$stnpair"].misc["dist"]

        if InputDict["method"] == "MWCS"
            time_axis, dt, err, coh = mwcs(ref,
                                           cur.corr,
                                           InputDict["freqmin"],
                                           InputDict["freqmax"],
                                           InputDict["fs"],
                                           InputDict["mintimelag"],
                                           InputDict["window_length"],
                                           InputDict["window_step"],
                                           InputDict["smoothing_half_win"])

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
            # set up time axis
            time_axis = InputDict["mintimelag"]:1/InputDict["fs"]:InputDict["maxtimelag"]

            window_length_samples = convert(Int,InputDict["window_length"] * InputDict["fs"])
            window_step_samples = convert(Int,InputDict["window_step"] * InputDict["fs"])
            minind = 1:window_step_samples:length(ref) - window_length_samples

            N = length(minind)
            dv_list = zeros(N)
            for ii=1:N
                window = collect(minind[ii]:minind[ii]+window_length_samples-1)
                dv, cc, cdp, eps, err, C = stretching(ref[:, 1],
                                                      cur.corr[:, 1],
                                                      time_axis,
                                                      window,
                                                      InputDict["freqmin"],
                                                      InputDict["freqmax"],
                                                      dvmin=-0.03,
                                                      dvmax=0.03)
                dv_list[ii] = dv
            end
            # average dv/V over all windows
            dvV=sum(dv_list)/N
        end
        # save dv/v to a GeoDataFrame (?)
    end
end

close(xcorrs)
close(ref)
