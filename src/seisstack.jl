include("utils.jl")
include("stacking.jl")
include("reference.jl")

using Statistics, Plots, FileIO, Printf, HashStack
export seisstack

"""
seisstack(InputDict::Dict)

compute stacking and save to jld2.
"""
function seisstack(InputDict::Dict)

	# #===
	# compute reference
	# ===#
	# # DEBUG
	#

	InputDict["Nfreqband"] = 10

	if InputDict["compute_reference"] == true
		if InputDict["stackmode"] == "linear" || InputDict["stackmode"] == "selective" ||  InputDict["stackmode"] == "hash"
			compute_reference_xcorr(InputDict)
			# error("selective reference has a bug in iteration. Currently not available.")
		elseif InputDict["stackmode"] == "robust"
			robust_reference_xcorr(InputDict)
		else
			#error("stackmode is either linear, selective or robust").
		end
	end

	#DEBUG:
	return 1
	#===
	compute stacking
	===#

    # load station pairs
    # compute stack in all station pairs by keys(ref)
	# if there is no reference, evaluate from InputDict["basefiname"]*".jld2"

	Output_rootdir = join(split(InputDict["basefiname"],"/")[1:end-3], "/") #.../OUTPUT
	Year = split(InputDict["basefiname"],"/")[end-2]

	reference = Output_rootdir*"/reference_xcorr_for$(Year).jld2" # this is fixed in the SeisXcorrelation/pstack.

	InputDict["referencepath"] = reference

    stapair = Array{String, 2}(undef, 3, 0)
    try
        f_ref = jldopen(reference, "r")
        pairs = split.(keys(f_ref), "-")
        pair1 = [pairs[x][1] for x = 1:length(pairs)]
		pair2 = [pairs[x][2] for x = 1:length(pairs)]
		comp = [pairs[x][3] for x = 1:length(pairs)]
        stapair = hcat(stapair, permutedims(hcat(pair1, pair2, comp)))
        close(f_ref)

    catch
        println("reference not found (if stackmode=linear, please ignore this warning.).\n
		 		use all potential station from basefile.")
	#	error("need to update. Abort")
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

    stapairlist = collect(zip(stapair[1,:], stapair[2,:], stapair[3,:]))
    InputDict["stapairlist"] = stapairlist

    # Parallelized by station pairs
    pmap(x -> map_stack(InputDict, x), stapairlist)

    str = "Stacking was successfully done.\n"
    printstyled(str; color=:green, bold=true)
	return nothing
end



"""
map_stack(InputDict, station)
stack process for pmap
"""
function map_stack(InputDict::Dict, station::Tuple)

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

    N_maxlag    = trunc(Int, 2*fs*maxlag + 1)

    if !haskey(InputDict, "filter")
        filter = false
    else
        filter = InputDict["filter"]
    end

	reference = InputDict["referencepath"]

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
	comp = station[3]

    stnpair = stn1*"-"*stn2*"-"*comp
    stnpairrev = stn2*"-"*stn1*"-"*comp

	# println(stnpair)
	# println(stnpairrev)

    #show progress
    progid = findfirst(x -> x==station, InputDict["stapairlist"])
    if mod(progid, 500) == 0
        println("station pairs No: $(progid)")
    end

    # output file name
    foname     = InputDict["fodir"]*"/stack_$(stn1)-$(stn2)-$(comp).jld2"

    # time lags for xcorr
    lags = -maxlag:1/fs:maxlag

	#save metadata
	xcorr_temp = get_metadata(timestamplist, timeslice, basefiname, N_maxlag, stnpair, stnpairrev)

	if typeof(xcorr_temp) != CorrData
		# no data available with this stnpair
		return nothing
	end

	# compute stack
    tscount = 0
	all_full_stnkeywithcha = String[]

    for (titer, time) in enumerate(timestamplist[timeslice[1]:timeslice[2]])
		if stackmode == "selective" || stackmode == "hash"
			# initialize removal fraction output
	        rmList = Array{Float64,1}()
	    end

        f_cur = jldopen(basefiname*".$time.jld2")
        full_stnkeys = try keys(f_cur[time]) catch; continue; end
		# reformat full_stnkeys to stack group name format
		nochan_stnkeys = String[]

		for k = 1:length(full_stnkeys)
			cur_stn1 = join(split(full_stnkeys[k], ".")[1:2], ".")
			cur_stn2 = join(split(full_stnkeys[k], ".")[5:6], ".")
			cur_comp = comp
			push!(nochan_stnkeys, cur_stn1*"-"*cur_stn2*"-"*cur_comp)
		end

		if stnpair ∈ nochan_stnkeys || stnpairrev ∈ nochan_stnkeys

            # # declare
            # xcorr = CorrData()
            # ref   = CorrData()

            if reference != false
				# read curent and reference
                if stnpair ∈ nochan_stnkeys
					ref_stnpair = stnpair
				else
					ref_stnpair = stnpairrev
				end

				fullstnpair_id = findfirst(x -> x==ref_stnpair, nochan_stnkeys)
				fullstnpair = full_stnkeys[fullstnpair_id]
                xcorr = try
					f_cur["$time/$fullstnpair"]
				catch
					println("current not found. make reference = false")
					reference = false
				end

				full_stnkeywithcha = fullstnpair

                f_ref = jldopen(reference, "r")

                ref = try
					f_ref[ref_stnpair]
				catch
					println("reference not found. make reference = false")
					reference = false

				end
                close(f_ref)

                if filter == true
                    xcorr.corr = bandpass(xcorr.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
                    ref.corr   = bandpass(ref.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
				end

            else
                # no reference:
				if stnpair ∈ nochan_stnkeys
					ref_stnpair = stnpair
				else
					ref_stnpair = stnpairrev
				end

				fullstnpair_id = findfirst(x -> x==ref_stnpair, nochan_stnkeys)
				fullstnpair = full_stnkeys[fullstnpair_id]
                xcorr = f_cur["$time/$fullstnpair"]
				full_stnkeywithcha = fullstnpair


                if filter == true
                    xcorr.corr = bandpass(xcorr.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
                end
            end

            nWins = length(xcorr.corr[1,:])

			xcorr.corr , nanzerocol = remove_nanandzerocol(xcorr.corr)
			xcorr.t = xcorr.t[nanzerocol]

			# stack shotttime-window cc per unit time

			if isempty(xcorr.corr) ||
				reference == false && stackmode == "selective" ||
				reference == false && stackmode == "hash"
				# skip this pair as there is no cc function
				# or no reference while selective/hash stack
				println("Current trace is empty. continue.")
				xcorr = CorrData()
				xcorr.fs = fs
				xcorr.maxlag = trunc(Int, N_maxlag)
				xcorr.corr = zeros(trunc(Int, N_maxlag),1)
				full_stnkeywithcha = ""

			# switch the stacking method
			elseif stackmode == "linear"

				stack!(xcorr, allstack=true)

            elseif stackmode == "selective"

				xcorr, ccList = selective_stack(xcorr, ref, InputDict)
				nRem = length(findall(x->(x<threshold), ccList))
				push!(rmList, nRem / nWins)

			elseif stackmode == "robust"

				robuststack!(xcorr)

			elseif stackmode == "hash"

				# @show xcorr.corr[200:400, 1]
				# @show ref.corr[200:400, 1]

				statsdir = InputDict["fodir"]*"/../hashstack_stats"
		        if !ispath(statsdir) mkpath(statsdir); end

				hashstack!(xcorr, ref,
					max_num_substack 		= 1e3,
					cc_substack_threshold 	= 0.3,
					SNR_threshold 			= 0.0,
					output_stats = true, # if output status hdf5 file to check the process
					output_stats_dir = statsdir, # directory of output status hdf5 file to check the process
					)

				cc_coef = cor(xcorr.corr, ref.corr)
				@show cc_coef

			else
				error("stack mode $stackmode does not exist.")

            end

            if stnpairrev ∈ nochan_stnkeys
                xcorr.corr = reverse(xcorr.corr, dims=1)
            end

            tscount += 1

        else
            # this station pair does not exist in this time stamp
            xcorr = CorrData()
            xcorr.fs = fs
            xcorr.maxlag = trunc(Int, N_maxlag)
            xcorr.corr = zeros(trunc(Int, N_maxlag),1)
			full_stnkeywithcha = ""
        end

        close(f_cur)

        norm_factor = maximum(abs.(xcorr.corr[:, 1]))

		if iszero(norm_factor)
			# this time is filled by zero; fill corr with NaN
			xcorr.corr .= NaN
			xcorr_temp.corr = hcat(xcorr_temp.corr, xcorr.corr[:, 1])

        elseif IsNormalizedampperUnit && !iszero(norm_factor) #to avoid deviding by zero
            xcorr_temp.corr = hcat(xcorr_temp.corr, (xcorr.corr[:, 1]./ norm_factor))
        else
            xcorr_temp.corr = hcat(xcorr_temp.corr, xcorr.corr[:, 1])
        end

		push!(all_full_stnkeywithcha, full_stnkeywithcha)

		# output selective removal fraction
		if (stackmode == "selective") && InputDict["IsOutputRemovalFrac"]
			#output removal fraction on this channel
			#println(rmList)
			if isempty(rmList)
				rmList = [NaN]
			end

			if stnpair ∈ nochan_stnkeys
				fname_out = join([time, stnpair, "selectiveremovalfraction","dat"], '.')
			else
				fname_out = join([time, stnpairrev, "selectiveremovalfraction","dat"], '.')
			end

			# y, jd = parse.(Int64, split(time, ".")[1:2])
			# tstamp_fname = replace(tstamp, ":" => "-")
			fodirtemp = split(basefiname, "/")
			fodir = join(fodirtemp[1:end-2], "/")*"/selectiveremoval_fraction"
			mkpath(fodir)
			fopath = fodir*"/"*fname_out
			open(fopath, "w") do io
			   write(io, @sprintf("%f\n", rmList[1]))
			end
		end
    end

	if tscount == 0
		#no data for this station pair
		#println("debug: nostationpair")
		return nothing
	end

    #===
    Manipulate stacking by unit num per stack
    ===#

    xcorr_all = deepcopy(xcorr_temp)
    xcorr_all.corr = Array{Float32, 2}(undef, trunc(Int, N_maxlag), 0)

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
		if Nstack == 0
			# This day ther is no stacked cc. fill NAN
			xcorr_all.corr = hcat(xcorr_all.corr, zeros(N_maxlag, 1) .* NaN)
        else
			xcorr_all.corr = hcat(xcorr_all.corr, xcorrstack ./ Nstack)
		end

        t1 = timestamplist[it]
        t2 = timestamplist[et]
        push!(T, get_midtime(t1, t2))
        #println(timestamp(get_midtime(t1, t2)))
        # slide with overlap
        icount += unitnumperstack - overlapperstack
    end

    xcorr_all.misc["T"] = T
	xcorr_all.misc["channelpair"] = all_full_stnkeywithcha

    # save xcorr
    jldopen(foname, "w") do file
        file["lags"] = lags
        file["xcorr"] = xcorr_all
        file["stackmode"] = stackmode
    end

    # dont plot if no cross-correlations were found. This would give div by 0
    if savefig
        # layout = Layout(width=1200, height=800,
        #                 xaxis_title="Lags [s]",
        #                 yaxis_title="Days since $(timestamplist[1])",
        #                 font=PlotlyJS.attr(size=11),
        #                 showlegend=false,
        #                 title="$stn1.$stn2")
		#
        # if tscount==0 println("No cross-correlations to plot. Skip."); end
        # trace = PlotlyJS.heatmap(x=lags, y=timestamp.(xcorr_all.misc["T"]), z=xcorr_all.corr, c=:RdBu)
        # p = PlotlyJS.Plot(trace, layout)
        # figdirtemp = split(basefiname, "/")
        # figdir = join(figdirtemp[1:end-2], "/")*"/fig"
        # if !ispath(figdir) mkpath(figdir); end
        # #figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)_$(timestamplist[1])."*figfmt
        # figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)."*figfmt
		# ORCA.restart_server()
		# ORCA.savefig(p, figname, format=figfmt)
		# println("test save data")
		if tscount==0 println("No cross-correlations to plot. Skip."); end
		Plots.pyplot()
		p = Plots.heatmap(lags, timestamp.(xcorr_all.misc["T"]),transpose(xcorr_all.corr), c=:balance)
		p = plot!(size=(800,600))
		figdirtemp = split(basefiname, "/")
		figdir = join(figdirtemp[1:end-2], "/")*"/fig"
		if !ispath(figdir) mkpath(figdir); end
		#figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)_$(timestamplist[1])."*figfmt
		figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)."*figfmt
		println(figname)
		Plots.savefig(p, figname)
		# save heatmap matrix for replot
		figjldname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2).jld2"
		FileIO.save(figjldname, "heatmap", xcorr_all.corr)
    end

    @show station
    return nothing
end

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
get_metadata(timestamplist::Array, timeslice::Array, basename::String, N_maxlag::Int)
get metadata from available cc jld2.
"""

function get_metadata(timestamplist::Array, timeslice::Array, basefiname::String, N_maxlag::Int,
	 	stnpair::String, stnpairrev::String)

	xcorr_temp = CorrData()

	breakflag=false
	for (titer, time) in enumerate(timestamplist[timeslice[1]:timeslice[2]])
		f_cur = jldopen(basefiname*".$time.jld2")
		full_stnkeys = keys(f_cur[time])
		# reformat full_stnkeys to stack group name format
		nochan_stnkeys = String[]

		for k = 1:length(full_stnkeys)
			cur_stn1 = join(split(full_stnkeys[k], ".")[1:2], ".")
			cur_stn2 = join(split(full_stnkeys[k], ".")[5:6], ".")
			cur_comp1 =split(full_stnkeys[k], ".")[4][end]
			cur_comp2 =split(full_stnkeys[k], ".")[8][end]
			cur_comp=cur_comp1*cur_comp2
			push!(nochan_stnkeys, cur_stn1*"-"*cur_stn2*"-"*cur_comp)
		end

		# println(stnpair)
		# println(stnpairrev)
		# println(nochan_stnkeys)
		#

		if stnpair ∈ nochan_stnkeys || stnpairrev ∈ nochan_stnkeys

			# declare: channel info will be lost through the time: only the first time step channel in CorrData is saved.

			xcorr = CorrData()

			if stnpair ∈ nochan_stnkeys
				ref_stnpair = stnpair
			else
				ref_stnpair = stnpairrev
			end

			if stnpair ∈ nochan_stnkeys
				fullstnpair_id = findfirst(x -> x==ref_stnpair, nochan_stnkeys)
				fullstnpair = full_stnkeys[fullstnpair_id]
				xcorr = f_cur["$time/$fullstnpair"]
			end
			# store meta data to xcorr_temp
			xcorr_temp = deepcopy(xcorr)
			xcorr_temp.corr = Array{Float32, 2}(undef, trunc(Int, N_maxlag), 0)

			return xcorr_temp
		end
	end

	# no data available with this stnpair
	return 1

end
