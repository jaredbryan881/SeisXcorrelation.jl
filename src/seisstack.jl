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

	timeunit   = InputDict["timeunit"]
	starttime  = InputDict["starttime"]
    endtime    = InputDict["endtime"]

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


    # station pairs
    stn1 = station[1]
    stn2 = station[2]
	comp = station[3]

    stnpair = stn1*"-"*stn2*"-"*comp
    stnpairrev = stn2*"-"*stn1*"-"*comp

    #show progress
    progid = findfirst(x -> x==station, InputDict["stapairlist"])
    if mod(progid, 500) == 0
        println("station pairs No: $(progid)")
    end

    # output file name
    foname = InputDict["fodir"]*"/stack_$(stn1)-$(stn2)-$(comp).jld2"

    # time lags for xcorr
    lags = -maxlag:1/fs:maxlag

	#save metadata
	xcorr_temp = get_metadata(timestamplist, starttime, endtime, basefiname, Nmaxlag, freqband,
							filter, stnpair, stnpairrev)

	if typeof(xcorr_temp) != CorrData
		# no data available with this stnpair
		return nothing
	end

	# intitialize time array
	xcorr_temp.t = Array{Float64, 1}(undef, 0)

	# compute stack
	all_full_stnkeywithcha = String[]

    for (titer, tstamp) in enumerate(timestamplist)

		y1, jd1 = parse.(Int64, split(tstamp, ".")[1:2])
	    m1, d1 = j2md(y1,jd1)
		datetime_tt = DateTime(y1, m1, d1)

		if  starttime > datetime_tt || endtime < datetime_tt
			continue;
		end

		if stackmode == "selective" || stackmode == "hash"
			# initialize removal fraction output
	        rmList = zeros(Float64, Nfreqband)
			rmList .= NaN
	    end

        f_cur = jldopen(basefiname*".$(tstamp).jld2")
        full_stnkeys = try keys(f_cur[tstamp]) catch; continue; end
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
					f_cur["$(tstamp)/$fullstnpair"]
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
                xcorr = f_cur["$(tstamp)/$fullstnpair"]
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

			append_wtcorr!(xcorr, freqband, α0 = InputDict["α0"], αmax = InputDict["αmax"])

			# stack shotttime-window cc per unit time
			for ifreq = 1:Nfreqband
				xcorr_ifreq = copy_without_wtcorr(xcorr)
				xcorr_ifreq.corr = xcorr.misc["wtcorr"][:,:,ifreq]

				if reference != false
					ref_ifreq = copy_without_wtcorr(ref)
					ref_ifreq.corr = reshape(ref.corr[:,ifreq], Nmaxlag, 1)
				end


				if isempty(xcorr_ifreq.corr) ||
					reference == false && stackmode == "selective" ||
					reference == false && stackmode == "hash"
					# skip this pair as there is no cc function
					# or no reference while selective/hash stack
					println("Current trace is empty or reference is not found. continue.")
					#xcorr_ifreq = CorrData()
					xcorr_ifreq.corr = zeros(trunc(Int, Nmaxlag),1)
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

					freqbandmin = xcorr.misc["freqband"][ifreq]
					freqbandmax = xcorr.misc["freqband"][ifreq+1]
					strfreq = @sprintf("%4.2f-%4.2f", round(freqbandmin, digits=2), round(freqbandmax, digits=2))

					statsdir = InputDict["Output_dir"]*"/hashstack_stats_"*strfreq
			        if !ispath(statsdir) mkpath(statsdir); end

					hashstack!(xcorr_ifreq, ref_ifreq,
						max_num_substack 		= 1e3,
						cc_substack_threshold 	= -1.0,
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

        else

			stacked_ifreq_cc .= NaN
			full_stnkeywithcha = ""
        end

        close(f_cur)

		# cat to all daily stack xcorr_temp

		wtcorr_reshaped = zeros(Float32, Nmaxlag, 1,Nfreqband)
		for ifreq = 1:Nfreqband
	        norm_factor = maximum(abs.(stacked_ifreq_cc[:, ifreq]), dims=1)[1]
			if iszero(norm_factor) || isnan(norm_factor)
				# this time is filled by zero; fill corr with NaN
				stacked_ifreq_cc[:, ifreq] .= NaN
				stacked_tr = stacked_ifreq_cc[:, ifreq]
	        elseif IsNormalizedampperUnit
				stacked_tr = stacked_ifreq_cc[:, ifreq] ./ norm_factor
	        else
				stacked_tr = stacked_ifreq_cc[:, ifreq]
	        end
			wtcorr_reshaped[:, 1, ifreq] = reshape(stacked_tr, length(stacked_ifreq_cc[:, ifreq]), 1, 1)
		end

		# append corr with freqband at this timestamp
		xcorr_temp.misc["wtcorr"] = cat(xcorr_temp.misc["wtcorr"],
											wtcorr_reshaped, dims=2)

		# append datetime in unixtime
		push!(xcorr_temp.t, d2u(datetime_tt))
		push!(all_full_stnkeywithcha, full_stnkeywithcha)

		# output selective removal fraction
		if (stackmode == "selective") && InputDict["IsOutputRemovalFrac"]
			#output removal fraction on this channel

			if stnpair ∈ nochan_stnkeys || stnpairrev ∈ nochan_stnkeys
				fname_out = join([tstamp, stnpair, "selectiveremovalfraction","dat"], '.')
			else
				continue;
			end

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



	if isempty(xcorr_temp.t)
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

	datetime_unittime = Dates.Second(timeunit)
	buffer_timewindow = 1.0 # [s]: shift time window backward to avoid including timestamp at right edge
	# initialize starttime
	twin_left  = starttime
	twin_right = twin_left + unitnumperstack*datetime_unittime

	icount = 0
	maxcount = 1e5

	while true
        #check if the window boundary is within start and end datetime.
		if twin_right > endtime
			# time window is outside of computed cc time
			break;
		end

		#find timestamp within time window
		datetime_list = u2d.(xcorr_temp.t)
		tind = findall(x -> twin_left-Second(buffer_timewindow) <= x &&
			 x < twin_right - Second(buffer_timewindow), datetime_list)

		if isempty(tind)
			# there is no time window within xcorr_temp. fill nan and slide to next time window.
			manip_nan_cc_reshaped = zeros(Nmaxlag, 1, Nfreqband) .* NaN
			xcorr_all.misc["wtcorr"] = cat(xcorr_all.misc["wtcorr"], manip_nan_cc_reshaped, dims=2)
			push!(T, (datetime2unix(twin_left) + datetime2unix(twin_right))/2)

			twin_left  = twin_right - overlapperstack*datetime_unittime
			twin_right = twin_left + unitnumperstack*datetime_unittime
			icount += 1
			if icount >= maxcount
				error("manipulate loop reaches limit $(maxcount). Please check overlapperstack, unitnumperstack, etc.")
			end
			continue;
		end

	 	manip_all_cc = zeros(Float32, Nmaxlag, Nfreqband)
		manip_temp_cc = @view xcorr_temp.misc["wtcorr"][:, tind, :]

		for ifreq = 1:Nfreqband
			manip_ifreq_cc = @view manip_temp_cc[:, :, ifreq]
			Nstack = 0
			# stack with each freqband
			manip_unitstacked_cc = zeros(Float32, Nmaxlag)
			for ii = 1:size(manip_ifreq_cc, 2)
				tr = manip_ifreq_cc[:,ii]
				if !any(isnan.(tr)) && !all(iszero.(tr))
					manip_unitstacked_cc .+= manip_ifreq_cc[:,ii]
					Nstack += 1
				end
			end

			if Nstack == 0
				manip_all_cc[:, ifreq] = zeros(Nmaxlag).*NaN
			else
				manip_all_cc[:, ifreq] = manip_unitstacked_cc ./ Nstack
			end
		end

		manip_all_cc_reshaped = reshape(manip_all_cc, Nmaxlag, 1, Nfreqband)
		xcorr_all.misc["wtcorr"] = cat(xcorr_all.misc["wtcorr"], manip_all_cc_reshaped, dims=2)
		push!(T, (datetime2unix(twin_left) + datetime2unix(twin_right))/2)

		# slide time window with overlap
		twin_left  = twin_right - overlapperstack*datetime_unittime
		twin_right = twin_left + unitnumperstack*datetime_unittime

		icount += 1
		if icount >= maxcount
			error("manipulate loop reaches limit $(maxcount). Please check overlapperstack, unitnumperstack, etc.")
		end
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
		if isempty(xcorr_temp.t); println("No cross-correlations to plot. Skip."); end
		Plots.pyplot()

		for ifreq = 1:Nfreqband
			freqbandmin = xcorr_all.misc["freqband"][ifreq]
			freqbandmax = xcorr_all.misc["freqband"][ifreq+1]

			xcorrplot = xcorr_all.misc["wtcorr"][:,:,ifreq]

			# normalize each timestamp by its maximum value
			for j = 1:size(xcorrplot, 2)
				norm_factor = maximum(abs.(xcorrplot[:, j]))[1]
				if iszero(norm_factor) || isnan(norm_factor)
					# this time is filled by zero; fill corr with NaN
					xcorrplot[:, j] .= NaN
				else
					xcorrplot[:, j] ./= norm_factor
				end
			end

			p = Plots.heatmap(lags, timestamp.(xcorr_all.misc["T"]),transpose(xcorrplot), c=:balance)
			p = Plots.plot!(size=(800,600),
					  title=@sprintf("%4.2f-%4.2f", round(freqbandmin, digits=2), round(freqbandmax, digits=2)),
					  xlabel = "Time lag [s]")

			figdir = InputDict["Output_dir"]*"/fig"
			if !ispath(figdir) mkpath(figdir); end
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

# """
# get_midtime(t1::String, t2::String)
# return average unix time between t1(ex. "2001.147.T00:00:00") and t2.
# """
# function get_midtime(t1::String, t2::String)
#     y1, jd1 = parse.(Int64, split(t1, ".")[1:2])
#     m1, d1 = j2md(y1,jd1)
#     time_init=DateTime(y1, m1, d1)
#
#     y2, jd2 = parse.(Int64, split(t2, ".")[1:2])
#     m2, d2 = j2md(y2,jd2)
#     time_end=DateTime(y2, m2, d2)
#     return (datetime2unix(time_init) + datetime2unix(time_end))/2
# end


"""
get_metadata(timestamplist::Array, timeslice::Array, basename::String, Nmaxlag::Int)
get metadata from available cc jld2.
"""

function get_metadata(timestamplist::Array, starttime::DateTime, endtime::DateTime, basefiname::String, Nmaxlag::Int,
	 	freqband::Union{Int, Array{Float64, 1}, Array{Float32, 1}}, filter::Union{Bool, AbstractArray},
		stnpair::String, stnpairrev::String)

	xcorr_temp = CorrData()

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

	for (titer, tstamp) in enumerate(timestamplist)

		y1, jd1 = parse.(Int64, split(tstamp, ".")[1:2])
		m1, d1 = j2md(y1,jd1)
		datetime_tt = DateTime(y1, m1, d1)

		if  starttime > datetime_tt || endtime < datetime_tt
			continue;
		end

		f_cur = jldopen(basefiname*".$(tstamp).jld2")
		full_stnkeys = keys(f_cur[tstamp])
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

		if stnpair ∈ nochan_stnkeys || stnpairrev ∈ nochan_stnkeys

			# declare: channel info will be lost through the time: only the first time step channel in CorrData is saved.

			xcorr = CorrData()

			if stnpair ∈ nochan_stnkeys
				ref_stnpair = stnpair
			else
				ref_stnpair = stnpairrev
			end

			fullstnpair_id = findfirst(x -> x==ref_stnpair, nochan_stnkeys)
			fullstnpair = full_stnkeys[fullstnpair_id]
			xcorr = f_cur["$(tstamp)/$fullstnpair"]

			# store meta data to xcorr_temp
			xcorr.corr , nanzerocol = remove_nanandzerocol(xcorr.corr)
			xcorr.t = xcorr.t[nanzerocol]

			if isempty(xcorr.corr) continue; end

			if filter != false
				xcorr.freqmin = filter[1]
				xcorr.freqmax = filter[2]
			end

			append_wtcorr!(xcorr, freqband)

			xcorr_temp = copy_without_wtcorr(xcorr)
			xcorr_temp.corr = Array{Float32, 2}(undef, Nmaxlag, 0)
			xcorr_temp.t = Array{Float64, 1}(undef, 0)
			xcorr_temp.misc["wtcorr"] = Array{Float32, 3}(undef, Nmaxlag, 0, Nfreqband)
			return xcorr_temp
		end
	end

	# no data available with this stnpair
	return 1

end
