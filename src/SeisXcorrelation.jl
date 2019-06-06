__precompile__()
module SeisXcorrelation

using SeisIO, SeisDownload, Noise, Printf, Dates, FFTW, JLD2, MPI

include("pairing.jl")
export seisxcorrelation


"""

    seisxcorrelation(maxtimelag::Real, finame::String, foname::String; IsAllComponent::Bool=false)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `finame::String,`    : Input file name e.g. network = "BPnetwork"
- `maxtimelag::Real`    : Maximum lag time e.g. maxtimelag = 100 [s]
- `corrtype::AbstractArray`    : Array of strings containing types of correlations to compute, e.g. ["xcorr", "acorr"]
- `corrorder`::Int    : Order of cross correlation, e.g. 3 for C3 (high order cross correlation)
- `IsAllComponent::Bool`   : If true, compute 3 X 3 components cross-correlation

# Output
- `foname.jld2`                 : contains SeisData structure with a hierarchical structure (CC function, metadata)

"""

function seisxcorrelation(finame::String, foname::String, corrtype::AbstractArray, corrorder::Int, maxtimelag::Real, freqmin::Real, freqmax::Real, fs::Real, cc_len::Int, cc_step::Int, IsAllComponents::Bool=false)

    MPI.Init()
    # establish the MPI communicator and obtain rank
    comm = MPI.COMM_WORLD
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    # read data from JLD2
    data = jldopen(finame)
    # output data file
    stlist = data["info/stationlist"]
    tstamplist = data["info/timestamplist"]

    # generate station pairs
    station_pairs = generate_pairs(stlist)
    # sort station pairs into autocorr, xcorr, and xchancorr
    sorted_pairs = sort_pairs(station_pairs)

    # create output file and save station information in JLD2
    if rank == 0
        jldopen(foname, "w") do file
            file["info/timestamplist"] = tstamplist;
            file["info/stationlist"] = stlist;
            file["info/sortedstationlist"] = sorted_pairs;
        end
    end

    mpiitrcount = 0

    baton = Array{Int32, 1}([0]) # for Relay Dumping algorithm

    for (t, tstamp) in enumerate(tstamplist)
        #parallelize one processor for one time stamp
        processID = t - (size * mpiitrcount)

        # if this mpiitrcount is final round or not
        length(stlist) - size >= size * mpiitrcount ? anchor_rank = size-1 : anchor_rank = mod(length(stlist), size)-1

        if rank == processID-1
            # iterate over correlation categories
            for ct in corrtype
                # iterate over station pairs in the current correlation category
                for j = 1:length(sorted_pairs[ct][1, :])
                    # get station names
                    stn1 = sorted_pairs[ct][1, j]
                    stn2 = sorted_pairs[ct][2, j]
                    println("Corrrelating $stn1 with $stn2")

                    # read station SeisChannels into SeisData before FFT
                    S1 = SeisData(data["$tstamp/$stn1"])
                    if ct=="acorr" S2=S1 else S2=SeisData(data["$tstamp/$stn2"]) end # S2 is a ref to S1 if "acorr"

                    # compute FFT using Noise.jl
                    FFT1 = compute_fft(S1, freqmin, freqmax, fs, cc_step, cc_len)
                    if ct=="acorr" FFT2=FFT1 else FFT2=compute_fft(S2, freqmin, freqmax, fs, cc_step, cc_len) end # FFT2 is a ref to FFT1 if "acorr"

                    # compute correlation using Noise.jl
                    xcorr = compute_cc(FFT1, FFT2, maxtimelag)

                    # save data after each cross correlation
                    varname = "$tstamp/$(xcorr.name)"
                    if size == 1
                        save_CorrData2JLD2(foname, varname, xcorr)
                        exit()

                    else
                        if rank == 0
                            save_CorrData2JLD2(foname, varname, xcorr)

                            if anchor_rank != 0
                                MPI.Send(baton, rank+1, 11, comm)
                                MPI.Recv!(baton, anchor_rank, 12, comm)
                            end

                        elseif rank == anchor_rank
                            MPI.Recv!(baton, rank-1, 11, comm)
                            save_CorrData2JLD2(foname, varname, xcorr)
                            MPI.Send(baton, 0, 12, comm)

                        else
                            MPI.Recv!(baton, rank-1, 11, comm)
                            save_CorrData2JLD2(foname, varname, xcorr)
                            MPI.Send(baton, rank+1, 11, comm)
                        end
                    end
                end
            end
        mpiitrcount += 1
        end
    end

    if rank == 0 println("Cross-correlation and data writing completed successfully.\nJob ended at "*string(now())) end

    MPI.Finalize()

    return 0
end

end

end
