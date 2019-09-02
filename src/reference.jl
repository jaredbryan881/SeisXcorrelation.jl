#include("utils.jl")
#include("io.jl")
include("stacking.jl")

export compute_reference_xcorr

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

    Output_rootdir = join(split(InputDict["basefiname"],"/")[1:end-3], "/") #.../OUTPUT

	refname = Output_rootdir*"/reference_xcorr.jld2" # this is fixed in the SeisXcorrelation/pstack.

	ref_dict_out = Dict()
	Numofrefstack = 0

	# get xcorr path between ref_starttime and ref_endtime
	ref_st = InputDict["ref_starttime"]
	ref_et = InputDict["ref_endtime"]
	ref_styear = Year(ref_st).value
	ref_etyear = Year(ref_et).value

	if InputDict["ref_iter"] < 1 error("ref_iter should be more than 1."); end

	for riter = 1:InputDict["ref_iter"] # iterate over "ref_iter" when using selective ref
		for year = ref_styear:ref_etyear
			corrname  = Output_rootdir*"/$(year)"*"/cc/$(year)_xcorrs"
			f = jldopen(corrname*".jld2");
			tslist = f["info/timestamplist"] # base xcorr file
			close(f)

			if riter == 1
				# first iteration should be stack="linear" because no reference
				ref_dicts = pmap(t->map_reference(t, InputDict, corr_name, stack="linear"), tslist)
				# collect all references into one dictionary at first iteration
				ref_dict_out_init = Dict()

				for i=1:length(ref_dicts)
					for stnpair in keys(ref_dicts[i])

						if haskey(ref_dict_out, stnpair)
							Numofrefstack += 1
							ref_dict_out_init[stnpair].corr .+= ref_dicts[i][stnpair].corr
						else
							# add new station pair into ref_dict_out
							ref_dict_out_init[stnpair] = copy(ref_dicts[i][stnpair])
						end
					end
				end

				# save initial reference
				f_out = jldopen(refname, "w")
				for stnpair in keys(ref_dict_out)
					f_out[stnpair] = ref_dict_out_init[stnpair]
				end
				close(f_out)

			else
				#allow either linear stack or selective stack
				ref_dicts = pmap(t->map_reference(t, InputDict, corr_name,
							stack=InputDict["stackmode"], reference=refname), tslist)
			end

			# merge all timestamp

			for i=1:length(ref_dicts)
				for stnpair in keys(ref_dicts[i])
					if haskey(ref_dict_out, stnpair)
						ref_dict_out[stnpair].corr .+= ref_dicts[i][stnpair].corr
						ref_dict_out[stnpair].misc["numofstack"] += 1

					else
						# initiate new station pair CorrData into ref_dict_out
						ref_dict_out[stnpair] = copy(ref_dicts[i][stnpair])
						ref_dict_out[stnpair].misc["numofstack"] = 1
					end
				end
			end

		end
	end

	# save final reference (this works even if riter = 1)
	f_out = jldopen(refname, "w")
	for stnpair in keys(ref_dict_out)
		f_out[stnpair] = ref_dict_out[stnpair]
	end
	close(f_out)

	# output reference status
	numofstackall = []
	for key = keys(ref_dict_out)
		push!(numofstackall, ref_dict_out[key].misc["numofstack"])
	end
	numofrefchannel = length(keys(ref_dict_out))
	maxnumofstack = maximum(numofstackall)
	minnumofstack = minimum(numofstackall)
	meannumofstack = mean(numofstackall)
	println("#---Reference stacking summary---#")
	println("Number of reference xcorr function: $(numofrefchannel)")
	println("Maximum num of stack: $(maxnumofstack)")
	println("Minimum num of stack: $(minnumofstack)")
	println("Mean num of stack   : $(meannumofstack)")
	println("#---reference xcorr is successfully saved---#\n
			$(refname)\n
			#--------------------------------------------#.")
end

"""
    map_reference(tstamp::String, InputDict::Dict, corrname::String;, stack::String="linear", reference::String="", thresh::Float64=-1.0)

Stack all cross-correlation functions for all given station pairs to generate a reference cross-correlation.

# Arguments
- `tstamp::String`    :
- `InputDict::Dict` : input dictionary
- `corrname::String,`    : Input base file name
- `phase_smoothing::Float64`    : Level of phase_smoothing (0 for linear stacking)
- `stack::String`     : "selective" for selective stacking and "linear" for linear stacking
- `reference::String`    : Path to the reference used in selective stacking
- `thresh::Float64`     : Threshold used for selective stacking

# Output
- `foname.jld2`    : contains arrays of reference cross-correlations for each station pair
"""
function map_reference(tstamp::String, InputDict::Dict, corrname::String; stackmode::String="linear", reference::String="")
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
	    grp = f_cur[tstamp] # xcorrs
	    println("$tstamp")

		# load reference
		if !isempty(reference)
			f_exref = jldopen(reference)
			ref = f_exref[pair]
			close(f_exref)
			if InputDict["filter"] !=false
				ref   = bandpass(ref.corr, InputDict["filter"][1], InputDict["filter"][2], xcorr.fs, corners=4, zerophase=false)
			end
		end

	    # iterate over station pairs
	    for pair in sort(keys(grp))
	        # TODO: use only unique station pairings when creating references. Currently no guarantee of uniqueness (reverse can exist)
	        # load xcorr
	        xcorr = try grp[pair] catch; continue end

	        #remove_nan!(xcorr)
	        # stack xcorrs over length of CorrData object using either "selective" stacking or "linear" stacking

			if InputDict["filter"] !=false
				xcorr = bandpass(xcorr.corr, InputDict["filter"][1], InputDict["filter"][2], xcorr.fs, corners=4, zerophase=false)
			end

	        if stackmode=="selective"
	            xcorr, rmList = selective_stacking(xcorr, ref, threshold=InputDict["threshold"],
								metric=InputDict["metric"], cohfilter = Inputdict["cohfilter"])
	        elseif stackmode=="linear"
				#linear stacking
	            stack!(xcorr, allstack=true, phase_smoothing=InputDict["phase_smoothing"])
			else
				error("stack mode $(stackmode) not defined.")
	        end

	        # stack xcorrs if they have a key, assign key if not
	        if haskey(ref_dict, pair)
	            ref_dict[pair].corr .+= xcorr.corr
	        else
	            ref_dict[pair] = copy(xcorr)
	        end
	    end

	    close(f_cur) # current xcorr file

	end

    return ref_dict
end
