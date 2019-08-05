using SeisIO, SeisNoise, JLD2, Statistics, FFTW, DSP, Plots

include("utils.jl")
include("partition.jl")

"""
    selective_stacking(data::CorrData, reference::CorrData; threshold::Float64=0.0, slice::Union{Bool, Float64, Array{Float64,1}}=false, metric::String="cc", win_len::Float64=10.0, win_step::Float64=5.0)

Stack the windows in a CorrData object that exceed a correlation-coefficient threshold with respect to a reference.

# Arguments
- `data::CorrData,`    : Input CorrData matrix used to define the reference and each window
- `reference::CorrData`    : Reference cross-correlation function to use in computing the correlation coefficient
- `threshold::Float64,`    : Value of correlation-coefficient below which we zero windows
- `slice::Union{Bool, Float64, Array{Float64,1}}`    : whether to slice the cross-correlations before computing convergence. If so, a single Float64 specifies the
                                                       lag at which to slice to keep just the ballistic wave and an array of Float64 specifies the start and end lag
                                                       to keep only the coda.

# Output
- `stackedData::CorrData,`    : Slectively stacked data
- `cList::Array{Float64,1}`    : Correlation coefficient of each window with respect to the reference
"""
function selective_stacking(data::CorrData, reference::CorrData; threshold::Float64=0.0, slice::Union{Bool, Float64, Array{Float64,1}}=false, metric::String="cc", coh_win_len::Float64=10.0, coh_win_step::Float64=5.0, filter::Union{Bool, Array{Float64,1}}=false)
    # slice data if given a time window
    # TODO: find a better way to pass arguments to this function that only apply to coh, or only to cc. e.g., win_len has no meaning if metric=="cc"
    if typeof(slice) != Bool
        ref = copy(reference)
        d = copy(data)
        if typeof(slice)==Float64
            # time vector
            tvec = -reference.maxlag:1/reference.fs:reference.maxlag
            # find all times within [-slice,slice] time
            t_inds = findall(x->(x>=-slice && x<=slice), tvec)

            # slice reference (maintain 2d array)
            ref.corr = reference.corr[t_inds, :]
            # slice data (maintain 2d array)
            d.corr = data.corr[t_inds, :]
        elseif typeof(slice)==Array{Float64,1}
            # convert startlag/endlag[s] to startlag/windowlength[samples]
            win_len = Int(diff(slice)[1] * reference.fs)
            startlag = Int(slice[1] * reference.fs)

            # partition reference to keep only the coda
            ref.corr = partition(reference.corr, startlag, win_len)
            # partition data to keep only the coda
            d.corr = partition(data.corr, startlag, win_len)
        else
            println("Please choose an allowable slicing operation. Exiting.")
            exit()
        end
    else
        # default to full reference cross-correlation
        ref = reference
        d = data
    end

    if metric=="cc"
        # compute correlation coefficient for each window in data with respect to reference
        cList = get_cc(d, ref)
    elseif metric=="coh"
        # compute coherence for each window as a function of frequency
        cohList = get_coh(d, ref, coh_win_len, coh_win_step)

        # use only coherencies for requested frequencies in the mean
        # default to using all frequencies
        if filter==false
            filter=[0.0, ref.fs/2]
        end
        # set acceptable frequencies to use in mean coh calculation
        freqs = rfftfreq(Int(coh_win_len*ref.fs), ref.fs)
        freq_inds = findall(x->(x>=filter[1] && x<=filter[2]), freqs)

        # compute correlation coefficient for each window in data with respect to reference
        ccList=get_cc(d, ref)

        # compute mean coherence
        mcohList = mean(cohList[freq_inds, :], dims=1)[1,:]

        # scale coherence by sign of the correlation coefficient to make extraction of constructively coherent traces possible
        cList = mcohList .* sign.(ccList)
    end

    if typeof(slice)==Array{Float64,1}
        # keep only xcorrs coherent in positive && negative coda
        good_fit_neg = findall(x->(x>=threshold), cList[1:Int(length(cList)/2)])
        good_fit_pos = findall(x->(x>=threshold), cList[Int(length(cList)/2+1):end])
        good_fit = intersect(good_fit_neg, good_fit_pos)
    else
        # find cross-correlations that fall below the cc threshold
        good_fit = findall(x->(x>=threshold), cList)
    end

    # copy data to extract well-fitting windows
    tempData = copy(data)
    tempData.corr = tempData.corr[:, good_fit]

    # linearly stack all data that exceeds the correlation-coefficient threshold
    stackedData = stack(tempData, allstack=true)

    plot=false
    if plot
        # coherence heat map
        p1=Plots.heatmap(1:length(cohList[1,:]), freqs, cohList)
        # Mean coherence and correlation coefficient
        p2 = Plots.plot(mcohList, label="Mean Coh", color= :red)
        plot!(p2, ccList, label="CC", color= :blue)
        plot!(p2, cList, label="Signed Coherence", color= :black, legend=false)

        # linear stack vs selective stack vs reference
        linstack = copy(data)
        linstack = stack(linstack, allstack=true)
        p3 = Plots.plot(normalize(linstack.corr[:,1]), label= "Linear")
        plot!(p3, normalize(stackedData.corr[:,1]), label="Selective")
        plot!(p3, normalize(ref.corr[:,1]), label="Reference")

        # spectrum of linear stack vs selectie stack vs reference
        fft_lin = rfft(linstack.corr[:,1])
        fft_sel = rfft(stackedData.corr[:,1])
        fft_ref = rfft(ref.corr[:, 1])
        freqs = rfftfreq(Int(length(ref.corr[:,1])*ref.fs), ref.fs)
        p4 = Plots.plot(freqs[1:500], normalize(abs.(real.(fft_lin)))[1:500], label="Linear")
        plot!(p4, freqs[1:500], normalize(abs.(real.(fft_sel)))[1:500], label="Selective")
        plot!(p4, freqs[1:500], normalize(abs.(real.(fft_ref)))[1:500], label="Reference")

        p5 = Plots.plot(p1, p2, layout=grid(2,1))
        p6 = Plots.plot(p3, p4, layout=grid(2,1))
        display(p5)
        readline()
        display(p6)
        readline()
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
        cList[j] = cor(reference.corr[:,1], data.corr[:,j])
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
