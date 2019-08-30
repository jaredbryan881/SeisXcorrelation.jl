include("./stacking.jl")

export pstack

using SeisIO, SeisNoise, JLD2, Dates, ORCA, PlotlyJS

"""
get_midtime(t1::String, t2::String)
return average unix time between t1(ex. "2001.147.T00:00:00") and t2.
"""
function get_midtime(t1::String, t2::String)
    y1, jd1 = parse.(Int64, split(t1, ".")[1:2])
    m1, d1 = j2md(y1,jd1)
    time_init=DateTime(y1, m1, d1)

    y2, jd2 = parse.(Int64, split(t2, ".")[1:2])
    m2, d2 = j2md(y2,jd2)
    time_end=DateTime(y2, m2, d2)
    return (datetime2unix(time_init) + datetime2unix(time_end))/2
end


"""
pstack(InputDict::Dict)

compute stacking and save to jld2.
"""
function pstack(InputDict::Dict)

    # load station pairs
    # If reference exists, suppose all station pairs are keys(ref)
    reference = InputDict["basefiname"]*"_ref.jld2"
    stapair = Array{String, 2}(undef, 2, 0)

    try
        f_ref = jldopen(reference, "r")
        pairs = split.(keys(f_ref), ".")
        pair1 = [join(pairs[x][1:4], ".") for x = 1:length(pairs)]
        pair2 = [join(pairs[x][5:8], ".") for x = 1:length(pairs)]
        stapair = hcat(stapair, permutedims(hcat(pair1, pair2)))
        close(f_ref)

    catch
        println("reference not found. use all potential station from basefile.")
        reference = false

        inFile = jldopen(InputDict["basefiname"]*".jld2", "r")
        corrsta         = inFile["info/corrstationlist"]
        close(inFile)

        corrtype   = InputDict["corrtype"]

        if any(occursin.("acorr", corrtype))
            stapair = hcat(stapair, corrsta["acorr"])
        end

        if any(occursin.("xcorr", corrtype))
            stapair = hcat(stapair, corrsta["xcorr"])
        end

        if any(occursin.("xcorrchan", corrtype))
            stapair = hcat(stapair, corrsta["xcorrchan"])
        end
    end

    datadir = split(InputDict["basefiname"], "/")
    fodir = join(datadir[1:end-2], "/")*"/stack"
    if ispath(fodir) rm(fodir, recursive=true);end
    mkpath(fodir)
    InputDict["fodir"] = fodir

    stapairlist = collect(zip(stapair[1,:], stapair[2,:]))
    InputDict["stapairlist"] = stapairlist

    #=== Debug ===#
    stapairlist = collect(zip(["NC.PADB..SHZ"], ["NC.PTA..SHZ"]))

    # Parallelized by station pairs
    pmap(x -> map_stack(InputDict, x), stapairlist)

    str = "Stacking was successfully done.\n"
    printstyled(str; color=:green, bold=true)
end



"""
map_stack(InputDict, station)
stack process for pmap
"""
function map_stack(InputDict::Dict, station::Tuple{String,String})

    basefiname = InputDict["basefiname"]
    stackmode  = InputDict["stackmode"] #"clustered_selective", "selective" or "linear"
    unitnumperstack = trunc(Int, InputDict["unitnumperstack"]) #num of stacking timestamp per unit stack (ex. 1 day = 1unit or 1month = 30unit)
    overlapperstack = trunc(Int, InputDict["overlapperstack"]) #num of overlapperstack per unit stack (ex. 1 day = 1unit or 1month = 30unit)
    savefig    = InputDict["savefig"]
    figfmt     = InputDict["figfmt"]
    IsNormalizedampperUnit   = InputDict["IsNormalizedampperUnit"]
    fs         = InputDict["fs"]
    maxlag     = InputDict["maxlag"]
    metric     = InputDict["metric"]
    threshold  = InputDict["threshold"]

    N_maxlag    = 2*fs*maxlag + 1

    if !haskey(InputDict, "filter")
        filter = false
    else
        filter = InputDict["filter"]
    end

    reference = basefiname*"_ref.jld2"
    if !ispath(reference)
        reference = false
    end

    inFile = jldopen(basefiname*".jld2", "r")
    timestamplist   = inFile["info/timestamplist"]
    close(inFile)

    # time index to stack
    if  InputDict["timeslice"] == false
        timeslice=[1,length(timestamplist)]
    else
        timeslice = InputDict["timeslice"]
    end

    # station pairs
    stn1 = station[1]
    stn2 = station[2]

    stnpair = stn1*"."*stn2
    stnpairrev = stn2*"."*stn1

    #show progress
    progid = findfirst(x -> x==station, InputDict["stapairlist"])
    if mod(progid, 500) == 0
        println("station pairs No: $(progid)")
    end

    # output file name
    foname     = InputDict["fodir"]*"/stack_$(stn1).$(stn2).jld2"

    # init variables
    xcorr_temp = CorrData()

    if stackmode == "clustered_selective"
        init_coef_series =  Array{Float64,1}()
        max_coef_series  =  Array{Float64,1}()
    elseif stackmode == "selective"
        rmList = Array{Float64,1}()
    end

    # time lags for xcorr
    lags = -maxlag:1/fs:maxlag

    tscount = 0
    for (titer, time) in enumerate(timestamplist[timeslice[1]:timeslice[2]])

        f_cur = jldopen(basefiname*".$time.jld2")
        stnkeys = keys(f_cur[time])

        # load cross-correlation or its reverse
        if stnpair ∈ stnkeys || stnpairrev ∈ stnkeys

            # declare
            xcorr = CorrData()
            ref   = CorrData()

            if reference!=false
                try
                    if stnpair ∈ stnkeys
                        xcorr = f_cur["$time/$stnpair"]
                        f_ref = jldopen(reference, "r")
                        ref = f_ref[stnpair]
                        close(f_ref)
                    else
                        xcorr = f_cur["$time/$stnpairrev"]
                        f_ref = jldopen(reference, "r")
                        ref = f_ref[stnpairrev]
                        close(f_ref)
                    end

                    if filter!=false
                        xcorr.corr = bandpass(xcorr.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
                        ref.corr   = bandpass(ref.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
                    end
                catch
                    println("reference not found. make reference = false")
                    reference = false
                end
            else
                # no reference:
                if stnpair ∈ stnkeys
                    xcorr = f_cur["$time/$stnpair"]
                else
                    xcorr = f_cur["$time/$stnpairrev"]
                end

                if filter!=false
                    xcorr.corr = bandpass(xcorr.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
                end
            end

            #if isempty(xcorr_temp.id)
	    	if xcorr_temp.fs == 0.0
                # store meta data to xcorr_temp
                xcorr_temp = deepcopy(xcorr)
                xcorr_temp.corr = Array{Float32, 2}(undef, trunc(Int, N_maxlag), 0)
            end

            nWins = length(xcorr.corr[1,:])

            if stackmode == "selective"

                if reference!=false
                    xcorr, ccList = selective_stacking(xcorr, ref, threshold=threshold, metric=metric, filter=filter)
                else
                    #xcorr, ccList = selective_stacking(xcorr, threshold=threshold, metric=metric, filter=filter)
                    stack!(xcorr, allstack=true, phase_smoothing=phase_smoothing)
                end

                nRem = length(findall(x->(x<threshold), ccList))
                push!(rmList, nRem / nWins)

            elseif stackmode == "clustered_selective"
                println("clustered_selective is under implementatin! do linear stacking")
                stack!(xcorr, allstack=true, phase_smoothing=phase_smoothing)

                # if reference!=false
                #     f_ref=jldopen(reference)
                #     ref = f_ref[stnpair]
                #
                #     if filter!=false
                #         ref.corr = bandpass(ref.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
                #     end
                #     #xcorr, ccList = selective_stacking(xcorr, ref, threshold=threshold, metric=metric, filter=filter, slice=slice)
                #     xcorr, initcoef, maxcoef = clustered_selective_stacking(xcorr, ref, threshold=threshold, metric=metric, filter=filter, slice=slice)
                #     close(f_ref)
                # end
                # #nRem = length(findall(x->(x<threshold), ccList))
                # #push!(rmList, nRem / nWins)
                # push!(init_coef_series, initcoef)
                # push!(max_coef_series, maxcoef)

            else
                stack!(xcorr, allstack=true, phase_smoothing=0.0)
            end

            if stnpairrev ∈ stnkeys
                xcorr.corr = reverse(xcorr.corr, dims=1)
            end

            tscount += 1

        else
            # this station pair does not exist in this time stamp

            xcorr = CorrData()
            xcorr.fs = fs
            xcorr.maxlag = trunc(Int, N_maxlag)
            xcorr.corr = zeros(trunc(Int, N_maxlag),1)

            #if isempty(xcorr_temp.id)
	    	if xcorr_temp.fs == 0.0
                # store meta data to xcorr_temp
                xcorr_temp = deepcopy(xcorr)
                xcorr_temp.corr = Array{Float32, 2}(undef, trunc(Int, N_maxlag), 0)
            end
        end

        close(f_cur)

        norm_factor = maximum(abs.(xcorr.corr[:, 1]))

        if IsNormalizedampperUnit && !iszero(norm_factor)
            xcorr_temp.corr = hcat(xcorr_temp.corr, (xcorr.corr[:, 1]./ norm_factor))
        else
            xcorr_temp.corr = hcat(xcorr_temp.corr, xcorr.corr[:, 1])
        end

    	if tscount == 0
	        #no data for this station pair
	        #println("debug: nostationpair")
			return nothing
		end
		
    end

    #===
    Manipulate stacking by unitnumperstack
    ===#

    xcorr_all = deepcopy(xcorr_temp)
    xcorr_all.corr = Array{Float32, 2}(undef, trunc(Int, N_maxlag), 0)

    corrdebug = xcorr_temp.corr

    if unitnumperstack <= overlapperstack
        error("unitnumperstack should be larger than overlapperstack.")
    end

    T = []
    icount=1
    while icount <= length(xcorr_temp.corr[1, :])-unitnumperstack+1

        it = icount
        et = icount + unitnumperstack - 1

        xcorrstack = zeros(N_maxlag)

        # stack over unitnumperstack
		Nstack = 0
        for ind = it:et
			if !any(isnan.(xcorr_temp.corr[:, ind])) && !all(iszero.(xcorr_temp.corr[:, ind]))
				xcorrstack .+= xcorr_temp.corr[:,ind]
				Nstack += 1

           	else
			   	#println("Nan or all zero found in corr. skip this unit time")
	   		end
        end
        #xcorr_all.corr = hcat(xcorr_all.corr, xcorrstack)
        xcorr_all.corr = hcat(xcorr_all.corr, xcorrstack ./ Nstack)

        t1 = timestamplist[it]
        t2 = timestamplist[et]
        push!(T, get_midtime(t1, t2))
        #println(timestamp(get_midtime(t1, t2)))
        # slide with overlap
        icount += unitnumperstack - overlapperstack
    end

    xcorr_all.misc["T"] = T

    # save xcorr
    jldopen(foname, "w") do file
        file["lags"] = lags
        file["xcorr"] = xcorr_all
        file["stackmode"] = stackmode
    end

    # dont plot if no cross-correlations were found. This would give div by 0
    if savefig
        layout = Layout(width=1200, height=800,
                        xaxis_title="Lags [s]",
                        yaxis_title="Days since $(timestamplist[1])",
                        font=PlotlyJS.attr(size=11),
                        showlegend=false,
                        title="$stn1.$stn2")

        if tscount==0 println("No cross-correlations to plot. Skip."); end
        trace = PlotlyJS.heatmap(x=lags, y=timestamp.(xcorr_all.misc["T"]), z=xcorr_all.corr, c=:RdBu)
        p = PlotlyJS.plot(trace, layout)
        figdirtemp = split(basefiname, "/")
        figdir = join(figdirtemp[1:end-2], "/")*"/fig"
        if !ispath(figdir) mkpath(figdir); end
        #figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)_$(timestamplist[1])."*figfmt
        figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)."*figfmt
        savefig(p, figname)
    end
end
