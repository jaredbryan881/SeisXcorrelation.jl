__precompile__()
module SeisXcorrelation

using SeisIO, Noise, Printf, Dates, FFTW, JLD2, MPI

export seisxcorrelation


"""

    seisxcorrelation(maxtimelag::Real, finame::String, foname::String; IsAllComponent::Bool=false)

Compute cross-correlation function and save data in jld2 file with SeisData format.

# Arguments
- `finame::String,`    : Input file name e.g. network = "BPnetwork"
- `maxtimelag::Real`    : Maximum lag time e.g. maxtimelag = 100 [s]
- `IsAllComponent::Bool`   : If true, compute 3 X 3 components cross-correlation

# Output
- `foname.jld2`                 : contains SeisData structure with a hierarchical structure (CC function, metadata)

"""
function seisxcorrelation(maxtimelag::Real, finame::String, foname::String; IsAllComponents::Bool=false) #Add whatever necessary arguments

    MPI.Init()
    # establish the MPI communicator and obtain rank
    comm = MPI.COMM_WORLD
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    #Read validation data from ["../EXAMPLE/codevalidation"]

    #Make Station pairs (after validating CC part, read from actual dataset and use statino list in SeisData structure)
    #Take into account IsAllComponents == true or not

    mpiitrcount = 0
    baton = Array{Int32, 1}([0]) # for Relay Dumping algorithm

    for i = 1:size(data["info/timestamplist"])

        #parallelize one processor for one time stamp
        processID = stid - (size * mpiitrcount)

        # if this mpiitrcount is final round or not
        length(stlist) - size >= size * mpiitrcount ? anchor_rank = size-1 : anchor_rank = mod(length(stlist), size)-1

        if rank == processID-1

            #roop each station pair

            for j = 1:length("""stationpairlist`""")

                # compute crosscorrelation using Noise.jl

                # store in SeisChannel structure
                # time, fs, misc["stationpair"], misc["cc1"]

                # save dataset

                # Save data one by one
                if size == 1
                    save_SeisData2JLD2(foname, varname, S)

                else
                    if rank == 0
                        save_SeisData2JLD2(foname, varname, S)

                        if anchor_rank != 0
                            MPI.Send(baton, rank+1, 11, comm)
                            MPI.Recv!(baton, anchor_rank, 12, comm)
                        end

                    elseif rank == anchor_rank
                        MPI.Recv!(baton, rank-1, 11, comm)
                        save_SeisData2JLD2(foname, varname, S)
                        MPI.Send(baton, 0, 12, comm)

                    else
                        MPI.Recv!(baton, rank-1, 11, comm)
                        save_SeisData2JLD2(foname, varname, S)
                        MPI.Send(baton, rank+1, 11, comm)
                    end
                end
            end

            mpiitrcount += 1

        end
    end

    if rank == 0 println("Crosscorrelation and Saving data is successfully done.\njob ended at "*string(now())) end

    MPI.Finalize()


    return 0
end

end
