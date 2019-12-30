using DSP

export selective_stack

"""
    selective_stack(data::CorrData, reference::CorrData, InputDict::Dict)

Stack the windows in a CorrData object that exceed a correlation-coefficient threshold with respect to a reference.

# Arguments
- `data::CorrData,`    : Input CorrData matrix used to define the reference and each window
- `reference::CorrData`    : Reference cross-correlation function to use in computing the correlation coefficient
- `InputDict::Dict,`    : Dictionary for input parameters


# Output
- `stackedData::CorrData,`    : Slectively stacked data
- `cList::Array{Float64,1}`    : Correlation coefficient of each window with respect to the reference
"""
function selective_stack(data::CorrData, reference::CorrData, InputDict::Dict)

	 threshold 			= InputDict["threshold"]
	 dtt_v 				= InputDict["dtt_v"]
	 metric 			= InputDict["metric"]
	 coh_win_len 		= InputDict["coh_win_len"]
	 coh_win_step 		= InputDict["coh_win_step"]
	 cohfilter 			= InputDict["cohfilter"]

    # slice data to remove ballistic wave if given ballistic wave velocity
    # TODO: find a better way to pass arguments to this function that only apply to coh, or only to cc. e.g., win_len has no meaning if metric=="cc"
    if typeof(dtt_v) != Bool
        ref = deepcopy(reference)
        d = deepcopy(data)
        # time vector

		ballistic_arrival_time = data.dist/dtt_v
		tvec = -data.maxlag:1/data.fs:data.maxlag

		if length(tvec) != size(d.corr, 1); error("length of time vector is not correct."); end

		coda_ind = findall(x->abs(x) >  ballistic_arrival_time, tvec)

        # slice reference (maintain 2d array)
        ref.corr = reference.corr[coda_ind, :]
        # slice data (maintain 2d array)
        d.corr = data.corr[coda_ind, :]

		# println(size(ref.corr, 1))
		# println(size(d.corr, 1))
		# @assert size(ref.corr, 1) == size(d.corr, 1)

    else
        # default to full reference cross-correlation
		ref = deepcopy(reference)
		d = deepcopy(data)
    end

    if metric=="cc"
        # compute correlation coefficient for each window in data with respect to reference
        cList = get_cc(d, ref)
    elseif metric=="coh"
        # compute coherence for each window as a function of frequency
        cohList = get_coh(d, ref, coh_win_len, coh_win_step)

        # use only coherencies for requested frequencies in the mean
        # default to using all frequencies
        if cohfilter==false
            cohfilter=[0.0, ref.fs/2]
        end
        # set acceptable frequencies to use in mean coh calculation
        freqs = DSP.rfftfreq(Int(coh_win_len*ref.fs), ref.fs)
        freq_inds = findall(x->(x>=cohfilter[1] && x<=cohfilter[2]), freqs)

        # compute correlation coefficient for each window in data with respect to reference
        ccList=get_cc(d, ref)

        # compute mean coherence
        mcohList = mean(cohList[freq_inds, :], dims=1)[1,:]

        # scale coherence by sign of the correlation coefficient to make extraction of constructively coherent traces possible
        cList = mcohList .* sign.(ccList)
    end

	good_fit = findall(x->(x>=threshold), cList)

    # copy data to extract well-fitting windows
	stackedData = deepcopy(data)
	stackedData.corr =	stackedData.corr[:, good_fit]

	# avoid NaN and zero column
	stackedData.corr , nanzerocol = remove_nanandzerocol(stackedData.corr)

	if isempty(stackedData.corr)
		# NaN padding as there is no reasonable stack
		stackedData.corr = zeros(length(data.corr[:,1]), 1)
		stackedData.corr .= NaN

	else
		# linearly stack all data that exceeds the correlation-coefficient threshold
		stack!(stackedData, allstack=true)
	end

    return stackedData, cList
end

"""
    get_cc(data::CorrData, reference::CorrData)

Compute the correlation coefficient for each window in a CorrData object with respect to a reference

# Arguments
- `data::CorrData`    : Input CorrData matrix used to define the reference and each window
- `reference::CorrData`    : Reference cross-correlation function to use in computing the correlation coefficient

# Output
- `cList::Array{Float64,1}`    :  Correlation coefficient for each window with respect to the reference
"""
function get_cc(data::CorrData, reference::CorrData)
    # array of correlation coefficients for each constituent cross-correlation window
    cList = Array{Float32,1}(undef,size(data.corr)[2])

    # iterate over each windowed cross-correlation and compute correlation-coefficient
    for j=1:size(data.corr)[2]
		try
        	cList[j] = cor(reference.corr[:,1], data.corr[:,j])
		catch y
			println(y)
			println(data)
			println(reference)
			cList[j] = -1.0
		end
    end

    return cList
end

"""
    get_coh(data::CorrData, reference::CorrData, win_len::Float64, win_step::Float64)

Compute the coherence between each window of a cross-correlation CorrData object with respect to a reference.

# Arguments
- `data::CorrData`    : Cross-correlation data to compute coherency
- `reference::CorrData`    : Reference cross-correlation function
- `win_len::Float64`    : Window length [s] when computing the psd using Welch's method. This parameter decides the resolvable frequencies in the coherence calculation
- `win_step::Float64`    : Step size [s] for windowing of each cross-correlation when computing the psd using Welch's method

# Outputs
- `coh::Array{Float32,2}`    : Matrix of size (nFreqs, nWindows) containing coherency values for each window and frequency
"""
function get_coh(data::CorrData, reference::CorrData, win_len::Float64, win_step::Float64)
    # convert window length and window step from time to samples
    win_len = Int(win_len * data.fs)
    win_step = Int(win_step * data.fs)

    # normalize signals
    data.corr = normalize(data.corr)
    reference.corr = normalize(reference.corr)

    # compute the coherence for each window in CorrData dim=2
    Ryy = welch_psd(reference.corr[:, 1], reference.corr[:, 1], data.fs, win_len, win_step)
    # initialize arrays to hold psd for each window
    Rxy = Array{eltype(data.corr), 2}(undef, size(Ryy,1), size(data.corr,2))
    Rxx = Array{eltype(data.corr), 2}(undef, size(Ryy,1), size(data.corr,2))
    # iterate over windows of CorrData
    for i=1:size(data.corr,2)
        # estimate the power spectral density using Welch's method
        Rxy[:, i] = welch_psd(data.corr[:, i], reference.corr[:, 1], data.fs, win_len, win_step, cross=true)
        Rxx[:, i] = welch_psd(data.corr[:, i], data.corr[:, i], data.fs, win_len, win_step)
    end

    # compute the sample coherence function
    coh = (Rxy.^2) ./ (Rxx .* Ryy)

    return coh
end

"""
    welch_psd(A::Array{R,1}, B::Array{R,1}, fs::Float64, win_len::Int64, win_step::Int64; cross=false) where R<:Real

Compute the cross power spectral density of two given arrays using Welch's method

Reference: Welch, P. D. (1967), "The use of Fast Fourier Transform for the estimation of power spectra:
A method based on time averaging over short, modified periodograms" (PDF),
IEEE Transactions on Audio and Electroacoustics, AU-15 (2): 70–73, doi:10.1109/TAU.1967.1161901

# Arguments
- `A::Array{R,1} where R<:Real`    : Signal 1
- `B::Array{R,1} where R<:Real`    : Signal 1
- `fs::Float64`    : sampling frequency
- `win_len::Int64`    : window length [samples] for FFT computation
- `win_step::Int64`    : step size [samples] for window advancement
- `cross::Bool`    : whether or not this is the psd of a cross-correlation

# Outputs
- `real.(cpsd)::Array{R,2}`    : Estimated cross power spectral density of two given time series
"""
function welch_psd(A::Array{R,1}, B::Array{R,1}, fs::Float64, win_len::Int64, win_step::Int64; cross=false) where R<:Real
    # partition the data into m windows
    A = window(A, win_len, win_step)
    B = window(B, win_len, win_step)

    # compute the FFT on the windowed data without bandpassing
    A_fft = process_fft(A, 0.0001, fs/2 - 0.0001, fs, corners=1)
    B_fft = process_fft(B, 0.0001, fs/2 - 0.0001, fs, corners=1)

    # compute the squared-magnitude of each signal block
    AB_fft2 = A_fft .* conj.(B_fft)

    # average over signal blocks to obtain an estimate of the cross power spectral density
    cpsd = mean(AB_fft2, dims=2)

    if cross
        power!(cpsd)
    end

    return real.(cpsd)
end

"""
    power!(A::AbstractArray)

Compute the power spectrum of a given FFT

# Arguments
- `A::AbstractArray`    : Input FFT

# Outputs
-  `power::AbstractArray`    : Power spectrum
"""
function power!(A::AbstractArray)
    return (abs.(A)).^2
end
power(A::AbstractArray) = (U = deepcopy(A);power!(U);return U)
