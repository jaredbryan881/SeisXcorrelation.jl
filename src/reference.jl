using Statistics

export compute_reference_xcorr, robust_reference_xcorr
"""
compute_reference_xcorr(InputDict::Dict)

compute and save reference cross-correlation function.
"""
function compute_reference_xcorr(InputDict::Dict)
    # computing reference xcorr
    #===
    Workflow:
    1. get timestamps between ref_starttime and ref_endtime and get path to xcorr files
    Note: only this process, you need to search among diferent years.
    2. At the first time, compute linear stack of all time stamp and sum up all of them
    3. iterate using selective stack
    4. save final reference dict
    ===#

	Input_rootdir = join(split(InputDict["basefiname"],"/")[1:end-3], "/") #.../OUTPUT

	Output_dir = InputDict["Output_dir"] #.../OUTPUT
	mkpath(Output_dir)

	refYear = split(InputDict["basefiname"],"/")[end-2]
	refname = Output_dir*"/../reference_xcorr_for$(refYear).jld2" # this is fixed in the SeisXcorrelation/pstack.

	ref_dict_out = Dict()

	# get xcorr path between ref_starttime and ref_endtime
	ref_st =InputDict["ref_starttime"]
	ref_et =InputDict["ref_endtime"]
	ref_styear = Dates.Year(ref_st).value
	ref_etyear = Dates.Year(ref_et).value


	# 1. make initial reference by linear stack
	ref_dict_dailystack = Dict()

	freqband = InputDict["freqband"]

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

	for year = ref_styear:ref_etyear

		corrname  = Input_rootdir*"/$(year)"*"/cc/$(year)_xcorrs"
		f = jldopen(corrname*".jld2");
		tslist = f["info/timestamplist"] # base xcorr file
		close(f)

		ref_dicts = []

		# first iteration should be stack="linear" because no reference
		ref_dicts = pmap(t->map_reference(t, InputDict, corrname, stackmode="linear"), tslist)
		#NOTE: ref_dicts has fillstation pass key including channels

		# Consider station pairs through the time regardless of channel
		# e.g. 1st day:NC.PDA..EHZ.NC.PCA..EHZ.  2nd day: NC.PCA..SHZ.NC.PDA..EHZ.
		# this case above we regard it as same station pair so that stacking togeter with "NC.PDA-NC.PCA-ZZ".

		for i=1:length(ref_dicts)
			for stnkey in keys(ref_dicts[i])
				xcorr_temp = ref_dicts[i][stnkey]
				stnkey = xcorr_temp.name
				stn1 = join(split(stnkey, ".")[1:2], ".")
				stn2 = join(split(stnkey, ".")[5:6], ".")
				comp = xcorr_temp.comp
				nochan_stnpair = stn1*"-"*stn2*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ
				nochan_stnpairrev = stn2*"-"*stn1*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ

				T = size(xcorr_temp.corr, 1)
				wtcorr_reshaped = reshape(xcorr_temp.corr, T, 1, Nfreqband)

				if haskey(ref_dict_dailystack, nochan_stnpair)
					if !isempty(xcorr_temp.corr)
	 					# ref_dict_dailystack[nochan_stnpair].corr = hcat(ref_dict_dailystack[nochan_stnpair].corr,
						#  ref_dicts[i][stnkey].corr)
						ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] =
						cat(ref_dict_dailystack[nochan_stnpair].misc["wtcorr"],
							wtcorr_reshaped, dims=2)
					end
		 		elseif haskey(ref_dict_dailystack, nochan_stnpairrev)
					if !isempty(xcorr_temp.corr)
						# ref_dict_dailystack[nochan_stnpair].corr = hcat(ref_dict_dailystack[nochan_stnpair].corr,
						#  reverse(ref_dicts[i][stnkey].corr, dims=1))
						ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] =
						cat(ref_dict_dailystack[nochan_stnpair].misc["wtcorr"],
							 reverse(wtcorr_reshaped, dims=1), dims=2)
					end
				else
					# add new station pair into ref_dict_out with stnpair (reversed stnpair in other time steps is taken into account above.)
					ref_dict_dailystack[nochan_stnpair] = deepcopy(xcorr_temp)
					ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] = wtcorr_reshaped
				end
			end
		end
	end

	# save initial reference by linear stack
	# DEBUG:
	@show refname
	f_out = jldopen(refname, "w")
	for stnpair in keys(ref_dict_dailystack)
		# #Debug
		#println(ref_dict_dailystack[stnpair])
		#@show any(isnan.(ref_dict_dailystack[stnpair].corr), dims=1)
		#@show any(isnan.(f_out[stnpair].corr), dims=1)
		# save initial reference by linear stack
		stacked_cc = copy_without_wtcorr(ref_dict_dailystack[stnpair])
		stacked_cc.corr = zeros(Float32, size(stacked_cc.corr, 1), Nfreqband)
		xcorr_ifreq = copy_without_wtcorr(ref_dict_dailystack[stnpair])

		for ifreq = 1:Nfreqband
			xcorr_ifreq.corr = ref_dict_dailystack[stnpair].misc["wtcorr"][:,:,ifreq]
			stack!(xcorr_ifreq, allstack=true)
			stacked_cc.corr[:, ifreq] = xcorr_ifreq.corr
		end
		# println(ref_dict_dailystack[stnpair])

		f_out[stnpair] = stacked_cc

	end
	close(f_out)


	if InputDict["stackmode"] == "linear" || InputDict["ref_iter"] == 0
		# end making reference
		return nothing
	end

	# iterate selective stack

	for riter = 1:InputDict["ref_iter"] # iterate over "ref_iter" when using selective ref
		println("reference method=selective iterate:$riter")

		ref_dict_dailystack = Dict()

		for year = ref_styear:ref_etyear

			corrname  = Input_rootdir*"/$(year)"*"/cc/$(year)_xcorrs"
			f = jldopen(corrname*".jld2");
			tslist = f["info/timestamplist"] # base xcorr file
			close(f)

			ref_dicts = []

			# first iteration should be stack="linear" because no reference
			#allow either linear stack or selective stack
			ref_dicts = pmap(t->map_reference(t, InputDict, corrname,
						stackmode="selective", reference=refname), tslist)

			for i=1:length(ref_dicts)
				for stnkey in keys(ref_dicts[i])
					xcorr_temp = ref_dicts[i][stnkey]
					stnkey = xcorr_temp.name
					stn1 = join(split(stnkey, ".")[1:2], ".")
					stn2 = join(split(stnkey, ".")[5:6], ".")
					comp = xcorr_temp.comp
					nochan_stnpair = stn1*"-"*stn2*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ
					nochan_stnpairrev = stn2*"-"*stn1*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ

					T = size(xcorr_temp.corr, 1)
					wtcorr_reshaped = reshape(xcorr_temp.corr, T, 1, Nfreqband)

					if haskey(ref_dict_dailystack, nochan_stnpair)
						if !isempty(xcorr_temp.corr)
		 					# ref_dict_dailystack[nochan_stnpair].corr = hcat(ref_dict_dailystack[nochan_stnpair].corr,
							#  ref_dicts[i][stnkey].corr)
							ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] =
							cat(ref_dict_dailystack[nochan_stnpair].misc["wtcorr"],
								wtcorr_reshaped, dims=2)
						end
			 		elseif haskey(ref_dict_dailystack, nochan_stnpairrev)
						if !isempty(xcorr_temp.corr)
							# ref_dict_dailystack[nochan_stnpair].corr = hcat(ref_dict_dailystack[nochan_stnpair].corr,
							#  reverse(ref_dicts[i][stnkey].corr, dims=1))
							ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] =
							cat(ref_dict_dailystack[nochan_stnpair].misc["wtcorr"],
								 reverse(wtcorr_reshaped, dims=1), dims=2)
						end
					else
						# add new station pair into ref_dict_out with stnpair (reversed stnpair in other time steps is taken into account above.)
						ref_dict_dailystack[nochan_stnpair] = deepcopy(xcorr_temp)
						ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] = wtcorr_reshaped
					end
				end
			end
		end

		# save reference by selective stack
		# update flowchart:
		# 1. change current reference name to be 'old'

		oldrefname = Output_dir*"/../old_reference_xcorr_for$(refYear).jld2" # this is fixed in the SeisXcorrelation/pstack.
		mv(refname, oldrefname, force=true)
		ref_in = jldopen(oldrefname, "r")
		f_out = jldopen(refname, "w")
		for stnpair in keys(ref_dict_dailystack)
			ref_old = ref_in[stnpair]

			stacked_cc = copy_without_wtcorr(ref_dict_dailystack[stnpair])
			stacked_cc.corr = zeros(Float32, size(stacked_cc.corr, 1), Nfreqband)
			xcorr_ifreq = copy_without_wtcorr(ref_dict_dailystack[stnpair])
			ref_ifreq = copy_without_wtcorr(ref_old)

			for ifreq = 1:Nfreqband

				xcorr_ifreq.corr = ref_dict_dailystack[stnpair].misc["wtcorr"][:,:,ifreq]
				ref_ifreq.corr = reshape(ref_old.corr[:,ifreq], length(ref_old.corr[:,ifreq]), 1)

				xcorr_ifreq = selective_stack(xcorr_ifreq, ref_ifreq, InputDict)[1]

				stacked_cc.corr[:, ifreq] = xcorr_ifreq.corr
			end
			# println(ref_dict_dailystack[stnpair])

			f_out[stnpair] = stacked_cc

		end
		close(ref_in)
		close(f_out)

	end

	println("#---selective stacking for reference xcorr is successfully saved---#\n$(refname)\n#--------------------------------------------#.")

	return nothing

end


"""
    map_reference(tstamp::String, InputDict::Dict, corrname::String; stackmode::String="linear", reference::String="")

Stack all cross-correlation functions for all given station pairs to generate a reference cross-correlation.

# Arguments
- `tstamp::String`    :
- `InputDict::Dict` : input dictionary
- `corrname::String,`    : Input base file name
- `stackmode::String`     : "selective" for selective stacking and "linear" for linear stacking
- `reference::String`    : Path to the reference used in selective stacking

"""
function map_reference(tstamp::String, InputDict::Dict, corrname::String; stackmode::String="linear", reference::String="")


	freqband = InputDict["freqband"]

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

	# hold reference xcorrs in memory and write all at once
	ref_dict = Dict()

	# if this time stamp is not in the ref_start and end time, skip stacking
	y, jd = parse.(Int64, split(tstamp, ".")[1:2])
	m, d  = j2md(y,jd)
	curdate=DateTime(y, m, d)

	if curdate >= InputDict["ref_starttime"] && curdate <= InputDict["ref_endtime"]
		#this time stamp is taken into account to stack

	    # read unstacked xcorrs for each time stamp
	    f_cur = try jldopen(corrname*".$tstamp.jld2")
			catch
				@show tstamp
				error("debug")
			end

	    grp = try
			f_cur[tstamp] # xcorrs
		catch
			ref_dict = Dict()
			return ref_dict
		end

	    println("$tstamp")

	    # iterate over station pairs
	    for pair in sort(keys(grp))
	        # Implemented -> TODO: use only unique station pairings when creating references. Currently no guarantee of uniqueness (reverse can exist)
	        # load xcorr
	        xcorr = try grp[pair] catch; continue end

			# remove nan and zero column
			xcorr.corr , nanzerocol = remove_nanandzerocol(xcorr.corr)
			xcorr.t = xcorr.t[nanzerocol]
			if isempty(xcorr.corr) continue; end

			if InputDict["filter"] !=false
				bandpass!(xcorr.corr, InputDict["filter"][1], InputDict["filter"][2], xcorr.fs, corners=4, zerophase=false)
				xcorr.freqmin = InputDict["filter"][1]
				xcorr.freqmax = InputDict["filter"][2]
			end

			# apply wavelet transform to C.corr and append it as 3D Array c.misc["wtcorr"]
			# basefiname = InputDict["basefiname"]
			# figdirtemp = split(basefiname, "/")
			#figdir = join(figdirtemp[1:end-2], "/")*"/fig_wtcorr"
			#figdir = InputDict["Output_dir"]
			figdir = ""
			append_wtcorr!(xcorr, freqband, figdir=figdir, α0 = InputDict["α0"], αmax = InputDict["αmax"])

			# load reference
			if stackmode=="selective"

				stnkey = xcorr.name
				stn1 = join(split(stnkey, ".")[1:2], ".")
				stn2 = join(split(stnkey, ".")[5:6], ".")
				comp = xcorr.comp
				nochan_stnpair = stn1*"-"*stn2*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ
				nochan_stnpairrev = stn2*"-"*stn1*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ

				f_exref = jldopen(reference)

				# NOTE: the try catch below keeps consitent to the order of station pairs between current and reference
				ref = try
						f_exref[nochan_stnpair]
					  catch
							try
								println("reversed station pair found.")
								f_exref[nochan_stnpairrev]
							catch
								println("debug: add new from second with linear.")
								stackmode="linear"
							end
					  end

				close(f_exref)

				if InputDict["filter"] !=false
					bandpass!(ref.corr, InputDict["filter"][1], InputDict["filter"][2], ref.fs, corners=4, zerophase=false)
					ref.freqmin = InputDict["filter"][1]
					ref.freqmax = InputDict["filter"][2]
				end

				if nochan_stnpairrev ∈ keys(f_exref)
					ref.corr = reverse(ref.corr, dims=1)
				end

			end

			# stack them and add to ref_dict.corr
			stacked_ifreq_cc = zeros(Float32, size(xcorr.corr, 1), Nfreqband)

			for ifreq = 1:Nfreqband
				xcorr_ifreq = copy_without_wtcorr(xcorr)
				xcorr_ifreq.corr = xcorr.misc["wtcorr"][:,:,ifreq]

				if stackmode=="linear"
					#linear stacking
					stack!(xcorr_ifreq, allstack=true)

				elseif stackmode=="selective"
					ref_ifreq = copy_without_wtcorr(ref)
					ref_ifreq.corr = reshape(ref.corr[:,ifreq], length(ref.corr[:,ifreq]), 1)

		            xcorr_ifreq, rmList = selective_stack(xcorr_ifreq, ref_ifreq, InputDict)

				else
					error("stack mode $(stackmode) not defined.")
		        end


				if ifreq == 1
					# initiate ref_dict metadata
					ref_dict[pair] = deepcopy(xcorr_ifreq)
				end

				# if there is no selected stack, skip this pair at this time
				if isempty(xcorr_ifreq.corr) || all(isnan.(xcorr_ifreq.corr)) continue; end
				stacked_ifreq_cc[:, ifreq] = xcorr_ifreq.corr

			end


	        # store all freqency bands in 2D CorrData.corr [:, Nfreqband]
			ref_dict[pair].corr = stacked_ifreq_cc
	    end

	    close(f_cur) # current xcorr file

	end

    return ref_dict
end

#=======================================================================#
#=======================================================================#
#=======================================================================#


"""
robust_reference_xcorr(InputDict::Dict)

compute and save reference cross-correlation function using robust stack.
"""
function robust_reference_xcorr(InputDict::Dict)
    # computing reference xcorr
    #===
    Workflow:
    1. get timestamps between ref_starttime and ref_endtime and get path to xcorr files
    Note: only this process, you need to search among diferent years.
    2. Compute day to day robust stack and store it into ref_corrdata
    3. Compute robust stack over reference period
    4. save final reference dict
    ===#

	Input_rootdir = join(split(InputDict["basefiname"],"/")[1:end-3], "/") #.../OUTPUT

	Output_dir = InputDict["Output_dir"] #.../OUTPUT
	if !ispath(Output_dir) mkpath(Output_dir); end

	refYear = split(InputDict["basefiname"],"/")[end-2]
	refname = Output_dir*"/../reference_xcorr_for$(refYear).jld2" # this is fixed in the SeisXcorrelation/pstack.

	ref_dict_out = Dict() # this contains {"stationpair" => CorrData}

	# get xcorr path between ref_starttime and ref_endtime
	ref_st =InputDict["ref_starttime"]
	ref_et =InputDict["ref_endtime"]
	ref_styear = Dates.Year(ref_st).value
	ref_etyear = Dates.Year(ref_et).value

	# collect all references into one dictionary at first iteration
	ref_dict_dailystack = Dict()

	freqband = InputDict["freqband"]

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

	for year = ref_styear:ref_etyear

		corrname  = Input_rootdir*"/$(year)"*"/cc/$(year)_xcorrs"
		f = jldopen(corrname*".jld2");
		tslist = f["info/timestamplist"] # base xcorr file
		close(f)

		ref_dicts = []

		# first iteration should be stack="linear" because no reference
		ref_dicts = pmap(t->map_robustreference(t, InputDict, corrname), tslist)
		#NOTE: ref_dicts has fillstation pass key including channels

		# Consider station pairs through the time regardless of channel
		# e.g. 1st day:NC.PDA..EHZ.NC.PCA..EHZ.  2nd day: NC.PCA..SHZ.NC.PDA..EHZ.
		# this case above we regard it as same station pair so that stacking togeter with "NC.PDA-NC.PCA-ZZ".

		for i=1:length(ref_dicts)
			for stnkey in keys(ref_dicts[i])
				xcorr_temp = ref_dicts[i][stnkey]
				stnkey = xcorr_temp.name
				stn1 = join(split(stnkey, ".")[1:2], ".")
				stn2 = join(split(stnkey, ".")[5:6], ".")
				comp = xcorr_temp.comp
				nochan_stnpair = stn1*"-"*stn2*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ
				nochan_stnpairrev = stn2*"-"*stn1*"-"*comp # e.g. NC.PDR-NC.PHA-ZZ

				T = size(xcorr_temp.corr, 1)
				wtcorr_reshaped = reshape(xcorr_temp.corr, T, 1, Nfreqband)

				if haskey(ref_dict_dailystack, nochan_stnpair)
					if !isempty(xcorr_temp.corr)
	 					# ref_dict_dailystack[nochan_stnpair].corr = hcat(ref_dict_dailystack[nochan_stnpair].corr,
						#  ref_dicts[i][stnkey].corr)
						ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] =
						cat(ref_dict_dailystack[nochan_stnpair].misc["wtcorr"],
							wtcorr_reshaped, dims=2)
					end
		 		elseif haskey(ref_dict_dailystack, nochan_stnpairrev)
					if !isempty(xcorr_temp.corr)
						# ref_dict_dailystack[nochan_stnpair].corr = hcat(ref_dict_dailystack[nochan_stnpair].corr,
						#  reverse(ref_dicts[i][stnkey].corr, dims=1))
						ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] =
						cat(ref_dict_dailystack[nochan_stnpair].misc["wtcorr"],
							 reverse(wtcorr_reshaped, dims=1), dims=2)
					end
				else
					# add new station pair into ref_dict_out with stnpair (reversed stnpair in other time steps is taken into account above.)
					ref_dict_dailystack[nochan_stnpair] = deepcopy(xcorr_temp)
					ref_dict_dailystack[nochan_stnpair].misc["wtcorr"] = wtcorr_reshaped
				end
			end
		end
	end


	# save final reference (this works even if riter = 1)
	f_out = jldopen(refname, "w")
	for stnpair in keys(ref_dict_dailystack)
		# #Debug
		#println(ref_dict_dailystack[stnpair])
		#@show any(isnan.(ref_dict_dailystack[stnpair].corr), dims=1)
		stacked_cc = copy_without_wtcorr(ref_dict_dailystack[stnpair])
		stacked_cc.corr = zeros(Float32, size(stacked_cc.corr, 1), Nfreqband)
		xcorr_ifreq = copy_without_wtcorr(ref_dict_dailystack[stnpair])


		for ifreq = 1:Nfreqband
			xcorr_ifreq.corr = ref_dict_dailystack[stnpair].misc["wtcorr"][:,:,ifreq]
			robuststack!(xcorr_ifreq)
			stacked_cc.corr[:, ifreq] = xcorr_ifreq.corr
		end
		# println(ref_dict_dailystack[stnpair])

		f_out[stnpair] = stacked_cc
		#@show any(isnan.(f_out[stnpair].corr), dims=1)

		# println(ref_dict_dailystack[stnpair])

	end
	close(f_out)

	# output reference status
	println("#---robust stacking for reference xcorr is successfully saved---#\n$(refname)\n#--------------------------------------------#.")
	return nothing
end

"""
    map_robustreference(tstamp::String, InputDict::Dict, corrname::String)

Robust stack cross-correlation functions for all given station pairs to generate a reference cross-correlation.

# Arguments
- `tstamp::String`    :
- `InputDict::Dict` : input dictionary
- `corrname::String,`    : Input base file name

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair
"""
function map_robustreference(tstamp::String, InputDict::Dict, corrname::String)

	freqband = InputDict["freqband"]

	if typeof(freqband) == Int
		Nfreqband = freqband
	else
		Nfreqband = length(freqband)-1
	end

    # hold reference xcorrs in memory and write all at once
	ref_dict = Dict()

	# if this time stamp is not in the ref_start and end time, skip stacking
	y, jd = parse.(Int64, split(tstamp, ".")[1:2])
	m, d  = j2md(y,jd)
	curdate=DateTime(y, m, d)

	if curdate >= InputDict["ref_starttime"] && curdate <= InputDict["ref_endtime"]
		#this time stamp is taken into account to stack

	    # read unstacked xcorrs for each time stamp
	    f_cur = jldopen(corrname*".$tstamp.jld2")
	    grp = try
			f_cur[tstamp] # xcorrs
		catch
			ref_dict = Dict()
			return ref_dict
		end

	    println("$tstamp")

	    # iterate over station pairs
	    for pair in sort(keys(grp))
	        # Implemented -> TODO: use only unique station pairings when creating references. Currently no guarantee of uniqueness (reverse can exist)
	        # load xcorr
	        xcorr = try grp[pair] catch; continue end

	        #remove_nan!(xcorr)
	        # stack xcorrs over length of CorrData object using either "selective" stacking or "linear" stacking

			xcorr.corr , nanzerocol = remove_nanandzerocol(xcorr.corr)
			xcorr.t = xcorr.t[nanzerocol]

			if isempty(xcorr.corr)
				#skip this pair as there is no cc function
				continue;
			end

			if InputDict["filter"] !=false
				bandpass(xcorr.corr, InputDict["filter"][1], InputDict["filter"][2], xcorr.fs, corners=4, zerophase=false)
				xcorr.freqmin = InputDict["filter"][1]
				xcorr.freqmax = InputDict["filter"][2]
			end

			# nancols = any(isnan.(xcorr.corr), dims=1)
			# xcorr.corr = xcorr.corr[:, vec(.!nancols)]

			# apply wavelet transform to C.corr and append it as 3D Array c.misc["wtcorr"]
			basefiname = InputDict["basefiname"]
			figdirtemp = split(basefiname, "/")
			#figdir = join(figdirtemp[1:end-2], "/")*"/fig_wtcorr"
			#figdir = InputDict["Output_dir"]
			figdir = ""

			append_wtcorr!(xcorr, freqband, figdir=figdir, α0 = InputDict["α0"], αmax = InputDict["αmax"])

			# @show any(isnan.(xcorr.corr), dims=1)
			# debugxcorr = deepcopy(xcorr)

			# println(vec(.!nancols))
			# println(xcorr.corr)

			# stack them and add to ref_dict.corr
			stacked_ifreq_cc = Array{Float32, 2}(undef, size(xcorr.corr, 1), Nfreqband)

			for ifreq = 1:Nfreqband
				xcorr_ifreq = copy_without_wtcorr(xcorr)
				xcorr_ifreq.corr = xcorr.misc["wtcorr"][:,:,ifreq]

				# TODO: implement slicing coda part
				#===
				# Draft code below. It will be made as function in util.jl

				min_ballistic_twin = InputDict["min_ballistic_twin"]
				dtt_v = InputDict["dtt_v"]
				max_coda_length = InputDict["max_coda_length"]
				tr = zeros(Float32, size(xcorr_ifreq.corr, 1))
				Ntimelag = length(tr/2)
				centerid = trunc(Int, Ntimelag/2)

				# minimum ballistic time window
	            min_ballistic_width = min_ballistic_twin * xcorr_ifreq.fs
	            minbal_window_left =  round(Int, centerid - min_ballistic_width)
	            minbal_window_right =  round(Int, centerid + min_ballistic_width)
	            minbal_window = vcat(collect(1:minbal_window_left), collect(window_right:Ntimelag))

	            if round(Int, xcorr_ifreq.fs * xcorr_ifreq.dist / dtt_v) < 1
	                # contains full ballistic part
	                coda_window_left =  round(Int, centerid - max_coda_length * xcorr_ifreq.fs)
	                coda_window_right =  round(Int, centerid + max_coda_length * xcorr_ifreq.fs)
	                window = intersect(collect(coda_window_left:coda_window_right),
	                                    minbal_window)
	            else
	                # remove ballistic wave part
	                coda_window_left =  round(Int, centerid - max_coda_length * xcorr_ifreq.fs) #[km] / [km/s]
	                coda_window_right =  round(Int, centerid + max_coda_length * xcorr_ifreq.fs) #[km] / [km/s]
	                ba_window_left = round(Int, centerid - xcorr_ifreq.fs * xcorr_ifreq.dist / dtt_v) #[km] / [km/s]
	                ba_window_right = round(Int, centerid + xcorr_ifreq.fs * xcorr_ifreq.dist / dtt_v) #[km] / [km/s]

	                if ba_window_left < 1
	                    ba_window_left = 1
	                end

	                if coda_window_left < 1
	                    coda_window_left = 1
	                end

	                if ba_window_right > Ntimelag
	                    ba_window_right = Ntimelag
	                end

	                if coda_window_right > Ntimelag
	                    coda_window_right = Ntimelag
	                end

	                window_left  = collect(coda_window_left:ba_window_left)
	                window_right = collect(ba_window_right:coda_window_right)
	                window = intersect(vcat(window_left, window_right),
	                                     minbal_window)
	                # check window
	                if minbal_window_left < minimum(window) || minbal_window_right > maximum(window)
	                    @warn("minimum ballistic window could be too large so that there is no stretching window.
	                        check balance of max_coda_length, attenuation_minthreshold and minimum ballistic window ")
	                end

	                if ba_window_left < coda_window_left || ba_window_right > coda_window_right || isempty(window)
	                    #this stationpare has no coda.
	                    continue;
	                end
	            end

				for it = 1:size(xcorr_ifreq.corr, 2)
					tr[window] = xcorr_ifreq.corr[window, it]
					xcorr_ifreq.corr[:, it] = tr
				end

				===#

				robuststack!(xcorr_ifreq)

				stacked_ifreq_cc[:, ifreq] = xcorr_ifreq.corr

				if ifreq == 1
					# initiate ref_dict metadata
					ref_dict[pair] = deepcopy(xcorr_ifreq)
				end

				# if there is no selected stack, skip this pair at this time
				if isempty(xcorr_ifreq.corr) || all(isnan.(xcorr_ifreq.corr)) continue; end

			end

			ref_dict[pair].corr = stacked_ifreq_cc

			# #
			# # print("nancheck:")
			# # @show any(isnan.(xcorr.corr), dims=1)
			#
			# if any(x-> x == true, any(isnan.(xcorr.corr), dims=1))
			# 	println("found nan in xcorr.corr.")
			# 	#@show(debugxcorr.corr)
			# 	robuststack_debug!(debugxcorr)
			#
			# 	#println(xcorr.corr)
			# 	error("Nan found in stacked trace. abort")
			# end

			# if there is no stacked trace, skip this pair

			# if isempty(xcorr.corr) || all(isnan.(xcorr.corr)) continue; end
			#
			# # stack xcorrs if they have a key, assign key if not
			# if haskey(ref_dict, pair)
			# 	ref_dict[pair].corr .+= xcorr.corr
			# else
			# 	ref_dict[pair] = deepcopy(xcorr)
			# end

	    end

	    close(f_cur) # current xcorr file

	end

    return ref_dict
end

#
#
# """
#   robuststack_debug(A)
# Performs robust stack on array `A`.
# Follows methods of Pavlis and Vernon, 2010.
# # Arguments
# - `A::AbstractArray`: Time series stored in columns.
# - `ϵ::AbstractFloat`: Threshold for convergence of robust stack.
# - `maxiter::Int`: Maximum number of iterations to converge to robust stack.
# """
# function robuststack_debug(A::AbstractArray{T};ϵ::AbstractFloat=Float32(1e-6),
#                      maxiter::Int=10) where T <: AbstractFloat
#     N = size(A,2)
#     Bold = median(A,dims=2)
#     w = Array{T}(undef,N)
#     r = Array{T}(undef,N)
#     d2 = Array{T}(undef,N)
#
#     # do 2-norm for all columns in A
#     for ii = 1:N
#         d2[ii] = norm(A[:,ii],2)
#     end
#
#     BdotD = sum(A .* Bold,dims=1)
#
#     for ii = 1:N
#
# 		@show BdotD[ii]
# 		@show Bold
# 		@show A[:,ii]
#
#         r[ii] = norm(A[:,ii] .- (BdotD[ii] .* Bold),2)
#         w[ii] = abs(BdotD[ii]) ./ d2[ii] ./ r[ii]
#
# 		@show r[ii]
# 		@show w[ii]
#
#     end
#
#     w ./= sum(w)
#
# 	@show r
#
# 	@show d2
#
#  	@show w
#
# 	@show A
#
#     Bnew = mean(A,weights(w),dims=2)
# 	@show Bnew
#
#     # check convergence
#     ϵN = norm(Bnew .- Bold,2) / (norm(Bnew,2) * N)
#     Bold = Bnew
#     iter = 0
#     while (ϵN > ϵ) && (iter <= maxiter)
#         BdotD = sum(A .* Bold,dims=1)
#
#         for ii = 1:N
#             r[ii] = norm(A[:,ii] .- (BdotD[ii] .* Bold),2)
#             w[ii] = abs(BdotD[ii]) ./ d2[ii] ./ r[ii]
#         end
#         w ./= sum(w)
#
#         Bnew = mean(A,weights(w),dims=2)
#
#         # check convergence
#         ϵN = norm(Bnew .- Bold,2) / (norm(Bnew,2) * N)
#         Bold = Bnew
#         iter += 1
#     end
#     return Bnew
# end
# robuststack_debug!(C::CorrData;ϵ::AbstractFloat=eps(Float32)) =
#        (C.corr = robuststack_debug(C.corr,ϵ=ϵ); C.t = C.t[1:1]; return C)
# robuststack_debug(C::CorrData;ϵ::AbstractFloat=eps(Float32))  =
#        (U = deepcopy(C); U.corr = robuststack_debug(U.corr,ϵ=ϵ); U.t = U.t[1:1];
#        return U)
