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
	#return 1
	#===
	compute stacking
	===#

    # load station pairs
    # compute stack in all station pairs by keys(ref)
	# if there is no reference, evaluate from InputDict["basefiname"]*".jld2"

	Input_rootdir = join(split(InputDict["basefiname"],"/")[1:end-3], "/") #.../OUTPUT

	Output_dir = InputDict["Output_dir"] #.../OUTPUT
	if !ispath(Output_dir) mkpath(Output_dir); end

	Year = split(InputDict["basefiname"],"/")[end-2]

	reference = Output_dir*"/../reference_xcorr_for$(Year).jld2" # this is fixed in the SeisXcorrelation/pstack.

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

    #datadir = split(InputDict["basefiname"], "/")
	#fodir = join(datadir[1:end-2], "/")*"/stack"
	fodir = Output_dir*"/stack"
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
	freqband   = InputDict["freqband"]

	Output_dir = InputDict["Output_dir"] #.../OUTPUT
	if !ispath(Output_dir) mkpath(Output_dir); end

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

    Nmaxlag    = trunc(Int, 2*fs*maxlag + 1)

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
	xcorr_temp = get_metadata(timestamplist, timeslice, basefiname, Nmaxlag, freqband,
							filter, stnpair, stnpairrev)

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
	        rmList = zeros(Float64, Nfreqband)
			rmList .= NaN
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

		stacked_ifreq_cc = zeros(Float32, Nmaxlag, Nfreqband)

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

				if filter != false && reference != false
					ref.corr    = bandpass(ref.corr, filter[1], filter[2], ref.fs, corners=4, zerophase=false)
					ref.freqmin = InputDict["filter"][1]
					ref.freqmax = InputDict["filter"][2]
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

            end

			nWins = length(xcorr.corr[1,:])

			xcorr.corr , nanzerocol = remove_nanandzerocol(xcorr.corr)
			xcorr.t = xcorr.t[nanzerocol]
			if isempty(xcorr.corr) continue; end
			if filter != false
				xcorr.corr    = bandpass(xcorr.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
				xcorr.freqmin = InputDict["filter"][1]
				xcorr.freqmax = InputDict["filter"][2]
			end

			append_wtcorr!(xcorr, freqband)

			# stack shotttime-window cc per unit time
			for ifreq = 1:Nfreqband
				xcorr_ifreq = copy_without_wtcorr(xcorr)
				xcorr_ifreq.corr = xcorr.misc["wtcorr"][:,:,ifreq]

				if reference != false
					ref_ifreq = copy_without_wtcorr(ref)
					ref_ifreq.corr = reshape(ref.corr[:,ifreq], Nmaxlag, 1)
				end


				if isempty(xcorr.corr) ||
					reference == false && stackmode == "selective" ||
					reference == false && stackmode == "hash"
					# skip this pair as there is no cc function
					# or no reference while selective/hash stack
					println("Current trace is empty. continue.")
					xcorr = CorrData()
					xcorr.fs = fs
					xcorr.maxlag = trunc(Int, Nmaxlag)
					xcorr.corr = zeros(trunc(Int, Nmaxlag),1)
					full_stnkeywithcha = ""

				# switch the stacking method
				elseif stackmode == "linear"

					stack!(xcorr_ifreq, allstack=true)

	            elseif stackmode == "selective"

					xcorr_ifreq, ccList = selective_stack(xcorr_ifreq, ref_ifreq, InputDict)
					nRem = length(findall(x->(x<threshold), ccList))
					rmList[ifreq] = nRem /  nWins

				elseif stackmode == "robust"

					robuststack!(xcorr_ifreq)

				elseif stackmode == "hash"

					# @show xcorr.corr[200:400, 1]
					# @show ref.corr[200:400, 1]

					statsdir = InputDict["Output_dir"]+"/hashstack_stats"
			        if !ispath(statsdir) mkpath(statsdir); end

					hashstack!(xcorr_ifreq, ref_ifreq,
						max_num_substack 		= 1e3,
						cc_substack_threshold 	= 0.3,
						SNR_threshold 			= 0.0,
						output_stats = true, # if output status hdf5 file to check the process
						output_stats_dir = statsdir, # directory of output status hdf5 file to check the process
						)

					cc_coef = cor(xcorr_ifreq.corr, ref_ifreq.corr)
					@show cc_coef

				else
					error("stack mode $stackmode does not exist.")
	            end

				if stnpairrev ∈ nochan_stnkeys
	                xcorr_ifreq.corr = reverse(xcorr_ifreq.corr, dims=1)
	            end
				# if there is no selected stack, skip this pair at this time
				if isempty(xcorr_ifreq.corr) || all(isnan.(xcorr_ifreq.corr)) continue; end

				stacked_ifreq_cc[:, ifreq] = xcorr_ifreq.corr

			end

			tscount += 1

        else
            # this station pair does not exist in this time stamp
            # xcorr = CorrData()
            # xcorr.fs = fs
            # xcorr.maxlag = trunc(Int, Nmaxlag)
            # xcorr.corr = zeros(trunc(Int, Nmaxlag),1)
			stacked_ifreq_cc .= NaN
			full_stnkeywithcha = ""
        end

        close(f_cur)

		# cat to all daily stack xcorr_temp

		wtcorr_reshaped = zeros(Float32, Nmaxlag, 1,Nfreqband)
		for ifreq = 1:Nfreqband
	        norm_factor = maximum(abs.(stacked_ifreq_cc[:, ifreq]), dims=1)
			if iszero(norm_factor)  #to avoid deviding by zero
				# this time is filled by zero; fill corr with NaN
				stacked_ifreq_cc[:, ifreq] .= NaN
				stacked_tr = stacked_ifreq_cc[:, ifreq]
				#xcorr_temp.corr = hcat(xcorr_temp.corr, xcorr.corr[:, 1])
	        elseif IsNormalizedampperUnit
	            #xcorr_temp.corr = hcat(xcorr_temp.corr, (xcorr.corr[:, 1]./ norm_factor))
				stacked_tr = stacked_ifreq_cc[:, ifreq] ./ norm_factor
	        else
	            #xcorr_temp.corr = hcat(xcorr_temp.corr, xcorr.corr[:, 1])
				stacked_tr = stacked_ifreq_cc[:, ifreq]
	        end
			wtcorr_reshaped[:, 1, ifreq] = reshape(stacked_tr, length(stacked_ifreq_cc[:, ifreq]), 1, 1)
		end
		xcorr_temp.misc["wtcorr"] = cat(xcorr_temp.misc["wtcorr"],
											wtcorr_reshaped, dims=2)

		push!(all_full_stnkeywithcha, full_stnkeywithcha)

		# output selective removal fraction
		if (stackmode == "selective") && InputDict["IsOutputRemovalFrac"]
			#output removal fraction on this channel
			#println(rmList)
			# if isempty(rmList)
			# 	rmList .= NaN
			# end

			if stnpair ∈ nochan_stnkeys || stnpairrev ∈ nochan_stnkeys
				fname_out = join([time, stnpair, "selectiveremovalfraction","dat"], '.')
			else
				continue;
			end

			# y, jd = parse.(Int64, split(time, ".")[1:2])
			# tstamp_fname = replace(tstamp, ":" => "-")
			# fodirtemp = split(basefiname, "/")
			# fodir = join(fodirtemp[1:end-2], "/")*"/selectiveremoval_fraction"
			fodir = InputDict["Output_dir"]*"/removal_fraction"
			mkpath(fodir)
			fopath = fodir*"/"*fname_out
			open(fopath, "w") do io
				for ii = 1:Nfreqband-1
					write(io, @sprintf("%f, ", rmList[ii]))
				end
				write(io, @sprintf("%f", rmList[Nfreqband]))
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
    xcorr_all.misc["wtcorr"] = Array{Float32, 3}(undef, trunc(Int, Nmaxlag), 0, Nfreqband)

    if unitnumperstack <= overlapperstack
        error("unitnumperstack should be larger than overlapperstack.")
    end

    T = []
    icount=1
    while icount <= length(xcorr_temp.misc["wtcorr"][1, :, 1])-unitnumperstack+1

        it = icount
        et = icount + unitnumperstack - 1

        # stack over unitnumperstack
		manup_ifreq_cc = zeros(Float32, Nmaxlag, Nfreqband)

		for ifreq = 1:Nfreqband
			xcorrstack = zeros(Nmaxlag)
			Nstack = 0
	        for ind = it:et
				if !any(isnan.(xcorr_temp.misc["wtcorr"][:, ind, ifreq])) && !all(iszero.(xcorr_temp.misc["wtcorr"][:, ind, ifreq]))
					xcorrstack .+= xcorr_temp.misc["wtcorr"][:, ind, ifreq]
					Nstack += 1
	           	else
				   	#println("Nan or all zero found in corr. skip this unit time")
		   		end
	        end

			if Nstack == 0
				manup_ifreq_cc[:, ifreq] = zeros(Nmaxlag, 1).*NaN
			else
				manup_ifreq_cc[:, ifreq] = xcorrstack ./ Nstack
			end
		end

		manupcorr_reshaped = reshape(manup_ifreq_cc, Nmaxlag, 1, Nfreqband)
		xcorr_all.misc["wtcorr"] = cat(xcorr_all.misc["wtcorr"], manupcorr_reshaped, dims=2)
        #xcorr_all.corr = hcat(xcorr_all.corr, xcorrstack)

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

		for ifreq = 1:Nfreqband
			freqbandmin = xcorr_all.misc["freqband"][ifreq]
			freqbandmax = xcorr_all.misc["freqband"][ifreq+1]

			xcorrplot = xcorr_all.misc["wtcorr"][:,:,ifreq]
			p = Plots.heatmap(lags, timestamp.(xcorr_all.misc["T"]),transpose(xcorrplot), c=:balance)
			p = Plots.plot!(size=(800,600),
					  title=@sprintf("%4.2f-%4.2f", round(freqbandmin, digits=2), round(freqbandmax, digits=2)),
					  xlabel = "Time lag [s]")

			# figdirtemp = split(basefiname, "/")
			# figdir = join(figdirtemp[1:end-2], "/")*"/fig"
			figdir = InputDict["Output_dir"]*"/fig"
			if !ispath(figdir) mkpath(figdir); end
			#figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)_$(timestamplist[1])."*figfmt
			strfreq = @sprintf("%4.2f-%4.2f", round(freqbandmin, digits=2), round(freqbandmax, digits=2))
			figname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)-"*strfreq*"."*figfmt
			println(figname)
			Plots.savefig(p, figname)
			# save heatmap matrix for replot
			figjldname = figdir*"/cc_$(stackmode)_$(stn1)-$(stn2)-"*strfreq*".jld2"
			FileIO.save(figjldname, Dict("heatmap" => xcorrplot, "T" => xcorr_all.misc["T"]))
		end
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
get_metadata(timestamplist::Array, timeslice::Array, basename::String, Nmaxlag::Int)
get metadata from available cc jld2.
"""

function get_metadata(timestamplist::Array, timeslice::Array, basefiname::String, Nmaxlag::Int,
	 	freqband::Union{Int, Array{Float64, 1}, Array{Float32, 1}}, filter::Union{Bool, AbstractArray},
		stnpair::String, stnpairrev::String)

	xcorr_temp = CorrData()

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

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
			xcorr.corr , nanzerocol = remove_nanandzerocol(xcorr.corr)
			xcorr.t = xcorr.t[nanzerocol]

			if isempty(xcorr.corr) continue; end
			if filter != false
				xcorr.corr    = bandpass(xcorr.corr, filter[1], filter[2], xcorr.fs, corners=4, zerophase=false)
				xcorr.freqmin = filter[1]
				xcorr.freqmax = filter[2]
			end

			append_wtcorr!(xcorr, freqband)

			xcorr_temp = copy_without_wtcorr(xcorr)
			xcorr_temp.corr = Array{Float32, 2}(undef, Nmaxlag, 0)
			xcorr_temp.misc["wtcorr"] = Array{Float32, 3}(undef, Nmaxlag, 0, Nfreqband)
			return xcorr_temp
		end
	end

	# no data available with this stnpair
	return 1

end
